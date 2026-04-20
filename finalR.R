# =============================================================================
# HOW DOES A TIME-WEIGHTED BRADLEY-TERRY MODEL, INTEGRATED WITH A BIVARIATE
# POISSON FRAMEWORK, COMPARE TO TRADITIONAL ELO RATINGS IN PREDICTING THE
# 2026 FIFA WORLD CUP CHAMPION?
# =============================================================================

# =============================================================================
# SECTION 1: Required Packages
# =============================================================================
library(worldfootballR)
library(dplyr)
library(BradleyTerry2)
library(ggplot2)
library(gt)
library(tidyr)
library(stringr)
library(lubridate)
library(MASS)
library(ggrepel)
library(viridis)
library(patchwork)
library(scales)
library(forcats)
library(tibble)

# =============================================================================
# SECTION 2: Raw Data
# =============================================================================
comps <- c(
  "FIFA World Cup",
  "FIFA World Cup Qualification — AFC",
  "FIFA World Cup Qualification — CAF",
  "FIFA World Cup Qualification — CONCACAF",
  "FIFA World Cup Qualification — CONMEBOL",
  "FIFA World Cup Qualification — UEFA",
  "UEFA Nations League",
  "UEFA European Football Championship",
  "UEFA European Football Championship Qualifying",
  "AFC Asian Cup",
  "Africa Cup of Nations",
  "CONCACAF Gold Cup",
  "CONMEBOL Copa América",
  "International Friendlies (M)"
)

df_raw <- load_match_comp_results(comp_name = comps)
colnames(df_raw)
print(head(df_raw$Home, 10))

# =============================================================================
# SECTION 3: DATA CLEANING AND PREPARATION
# =============================================================================

# Team name cleaning
clean_team_name <- function(name) {
  name <- str_trim(name)
  # Remove leading lowercase 2-letter code + space: "us United States"
  name <- str_replace(name, "^[a-z]{2}\\s+", "")
  # Remove trailing space + lowercase 2-letter code: "United States us"
  name <- str_replace(name, "\\s+[a-z]{2}$", "")
  # Safety: also strip any UPPERCASE 2-letter code variants
  name <- str_replace(name, "^[A-Z]{2}\\s+", "")
  name <- str_replace(name, "\\s+[A-Z]{2}$", "")
  str_squish(name)
}

country_map <- c(
  "USA"                  = "United States",
  "Korea Republic"       = "South Korea",
  "Korea Rep"            = "South Korea",
  "Republic of Korea"    = "South Korea",
  "Iran"                 = "IR Iran",
  "Ivory Coast"          = "Cote d'Ivoire",
  "Cape Verde"           = "Cabo Verde",
  "Czech Republic"       = "Czechia",
  "Bosnia"               = "Bosnia and Herzegovina",
  "Trinidad"             = "Trinidad and Tobago",
  "United Arab Emirates" = "UAE",
  "China PR"             = "China",
  "Chinese Taipei"       = "Taiwan",
  "DPR Korea"            = "North Korea",
  "eng England"          = "England",
  "england"              = "England",
  "England"              = "England"
)

standardise_name <- function(name) {
  cleaned <- clean_team_name(name)
  if (cleaned %in% names(country_map)) return(unname(country_map[cleaned]))
  cleaned
}

# Other factors
df_clean <- df_raw %>%
  rename(Competition = Competition_Name) %>%
  filter(!is.na(HomeGoals), !is.na(AwayGoals),
         !is.na(Home),      !is.na(Away)) %>%
  mutate(
    Date   = as.Date(Date),
    team1  = sapply(Home, standardise_name),
    team2  = sapply(Away, standardise_name),
    score1 = as.integer(HomeGoals),
    score2 = as.integer(AwayGoals)
  ) %>%
  filter(Date >= as.Date("2018-01-01")) %>%     # post-2018 window
  filter(nchar(team1) > 2, nchar(team2) > 2) %>% # drop unresolved codes
  distinct(Date, team1, team2, score1, score2, .keep_all = TRUE) %>%
  mutate(
    total_goals = score1 + score2,
    goal_diff   = score1 - score2,
    result      = case_when(
      score1 > score2 ~ "home_win",
      score1 < score2 ~ "away_win",
      TRUE            ~ "draw"
    )
  )

cat(sprintf("After initial cleaning: %d matches\n", nrow(df_clean)))
print(head(sort(unique(c(df_clean$team1, df_clean$team2))), 20))

# Competition-importance weights
comp_weights <- c(
  "FIFA World Cup"                                 = 3.0,
  "CONMEBOL Copa America"                          = 2.5,
  "UEFA European Football Championship"            = 2.5,
  "Africa Cup of Nations"                          = 2.0,
  "AFC Asian Cup"                                  = 2.0,
  "CONCACAF Gold Cup"                              = 2.0,
  "FIFA World Cup Qualification - UEFA"            = 1.8,
  "FIFA World Cup Qualification - CONMEBOL"        = 1.8,
  "FIFA World Cup Qualification - AFC"             = 1.5,
  "FIFA World Cup Qualification - CAF"             = 1.5,
  "FIFA World Cup Qualification - CONCACAF"        = 1.5,
  "UEFA Nations League"                            = 1.5,
  "UEFA European Football Championship Qualifying" = 1.3,
  "International Friendlies (M)"                   = 0.5
)

df_clean <- df_clean %>%
  mutate(
    comp_weight = ifelse(Competition %in% names(comp_weights),
                         comp_weights[Competition], 1.0)
  )

# Time-decay weights
# Exponential decay, half-life ~2.5 years (912 days)
t_now        <- max(df_clean$Date, na.rm = TRUE)
lambda_decay <- log(2) / 912
df_clean <- df_clean %>%
  mutate(
    days_ago     = as.numeric(t_now - Date),
    time_weight  = exp(-lambda_decay * days_ago),
    final_weight = time_weight * comp_weight
  ) %>%
  mutate(final_weight = final_weight / max(final_weight))  # normalise to [0,1]


# Sparsity filter: keep teams with >= 8 total matches
team_match_counts <- df_clean %>%
  pivot_longer(c(team1, team2), values_to = "team") %>%
  count(team) %>%
  filter(n >= 8)
eligible_teams <- sort(team_match_counts$team)

df_clean <- df_clean %>%
  filter(team1 %in% eligible_teams, team2 %in% eligible_teams)
all_teams <- sort(unique(c(df_clean$team1, df_clean$team2)))

# Bradley-Terry subset (decisive matches only)
df_bt_all <- df_clean %>% filter(result != "draw")
bt_counts <- df_bt_all %>%
  pivot_longer(c(team1, team2), values_to = "team") %>%
  count(team) %>%
  filter(n >= 5)

bt_teams <- sort(bt_counts$team)
df_bt_filtered <- df_bt_all %>%
  filter(team1 %in% bt_teams, team2 %in% bt_teams) %>%
  mutate(
    team1   = factor(team1, levels = bt_teams),
    team2   = factor(team2, levels = bt_teams),
    bt_win1 = as.integer(score1 > score2),  # 1 if team1 won, else 0
    bt_win2 = as.integer(score2 > score1)   # 1 if team2 won, else 0
  )

cat(sprintf("BT dataset: %d decisive matches | %d teams\n",
            nrow(df_bt_filtered), length(bt_teams)))

# =============================================================================
# SECTION 4: TIME-WEIGHTED BRADLEY-TERRY MODEL
# =============================================================================
levels(df_bt_filtered$team1)[grep("(?i)england", levels(df_bt_filtered$team1), perl = TRUE)] <- "England"
levels(df_bt_filtered$team2)[grep("(?i)england", levels(df_bt_filtered$team2), perl = TRUE)] <- "England"

bt_model <- BTm(
  outcome  = cbind(bt_win1, bt_win2),
  player1  = team1,
  player2  = team2,
  weights  = final_weight,
  data     = df_bt_filtered,
  id       = "team"
)

raw_abilities <- BTabilities(bt_model)
bt_abilities <- as.data.frame(raw_abilities) %>%
  rownames_to_column("team") %>%
  rename(se_bt = s.e.) %>%
  arrange(desc(ability)) %>%
  mutate(
    rank_bt      = row_number(),
    win_prob_avg = plogis(ability)  # P(beat an average opponent)
  )

# Use 2 decimal places
bt_abilities_formatted <- bt_abilities %>%
  mutate(
    ability = round(ability, 2),
    se_bt = round(se_bt, 2),
    win_prob_avg = round(win_prob_avg, 2)
  )

#Results & Save
print(head(bt_abilities_formatted[, c("rank_bt", "team", "ability", "se_bt")], 20))
print(colnames(bt_abilities_formatted))

top20_table <- bt_abilities_formatted %>%
  slice_head(n = 20) %>%
  dplyr::select(rank_bt, team, ability, se_bt) %>%
  gt() %>%
  tab_header(
    title = "Top 20 Teams by Bradley-Terry Ability",
    subtitle = "Time-weighted estimates with standard errors (2018-present)"
  ) %>%
  cols_label(
    rank_bt = "Rank",
    team = "Team",
    ability = "BT Ability",
    se_bt = "Std. Error"
  ) %>%
  fmt_number(
    columns = c(ability, se_bt),
    decimals = 2
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = rank_bt, rows = rank_bt <= 3)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgoldenrodyellow"),
    locations = cells_body(rows = rank_bt <= 3)
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 12,
    table.width = pct(80)
  )

gtsave(top20_table, "C:\\Users\\freya\\Downloads\\SMGT430\\top20_bt_table.pdf")

# =============================================================================
# SECTION 5: BIVARIATE POISSON FRAMEWORK
# =============================================================================

# Prepare data - fix England name first
df_clean$team1[df_clean$team1 == "eng England"] <- "England"
df_clean$team2[df_clean$team2 == "eng England"] <- "England"
df_clean$team1[df_clean$team1 == "england"] <- "England"
df_clean$team2[df_clean$team2 == "england"] <- "England"

df_pois <- df_clean %>%
  filter(team1 %in% eligible_teams, team2 %in% eligible_teams) %>%
  mutate(
    score1 = pmin(score1, 8L),
    score2 = pmin(score2, 8L)
  )

# Calculate simple goal-based ratings
pois_teams <- sort(unique(c(df_pois$team1, df_pois$team2)))
n_teams <- length(pois_teams)

# Calculate goal difference per match for each team
team_goals <- data.frame()
for(team in pois_teams) {
  # Matches as home
  home_matches <- df_pois[df_pois$team1 == team, ]
  # Matches as away  
  away_matches <- df_pois[df_pois$team2 == team, ]
  
  goals_for <- sum(home_matches$score1) + sum(away_matches$score2)
  goals_against <- sum(home_matches$score2) + sum(away_matches$score1)
  matches <- nrow(home_matches) + nrow(away_matches)
  
  team_goals <- rbind(team_goals, data.frame(
    team = team,
    goals_for = goals_for,
    goals_against = goals_against,
    matches = matches,
    stringsAsFactors = FALSE
  ))
}

# Calculate simple ratings
team_goals$goal_diff <- team_goals$goals_for - team_goals$goals_against
team_goals$goal_diff_per_match <- team_goals$goal_diff / team_goals$matches
# Scale to z-scores
team_goals$net_strength <- scale(team_goals$goal_diff_per_match)[,1]
team_goals$attack_strength <- team_goals$net_strength / 2
team_goals$defence_strength <- -team_goals$net_strength / 2
# Sort and add rank
team_goals <- team_goals[order(team_goals$net_strength, decreasing = TRUE), ]
team_goals$rank_pois <- 1:nrow(team_goals)
# Round numeric columns to 2 decimal places
team_goals$attack_strength <- round(team_goals$attack_strength, 2)
team_goals$defence_strength <- round(team_goals$defence_strength, 2)
team_goals$net_strength <- round(team_goals$net_strength, 2)

# Create pois_abilities dataframe
pois_abilities <- team_goals[, c("rank_pois", "team", "attack_strength", 
                                 "defence_strength", "net_strength")]
pois_abilities$team[pois_abilities$team == "England eng"] <- "England"

# Calculate a simple NLL approximation
calculate_simple_nll <- function() {
  total_ll <- 0
  for(i in 1:nrow(df_pois)) {
    team1_idx <- which(pois_abilities$team == df_pois$team1[i])
    team2_idx <- which(pois_abilities$team == df_pois$team2[i])
    if(length(team1_idx) == 0 || length(team2_idx) == 0) next
    team1_rating <- pois_abilities$net_strength[team1_idx]
    team2_rating <- pois_abilities$net_strength[team2_idx]
    
    # Expected goals based on ratings
    exp_goals1 <- exp(0.5 + team1_rating - team2_rating)
    exp_goals2 <- exp(0.5 + team2_rating - team1_rating)
    
    # Poisson log-likelihood
    ll1 <- dpois(df_pois$score1[i], exp_goals1, log = TRUE)
    ll2 <- dpois(df_pois$score2[i], exp_goals2, log = TRUE)
    total_ll <- total_ll + (ll1 + ll2) * df_pois$final_weight[i]
  }
  return(-total_ll / nrow(df_pois))
}

nll_approx <- calculate_simple_nll()

#Results & save
cat(sprintf("Approximate Negative Log-Likelihood: %.2f\n", nll_approx))
cat(sprintf("Teams evaluated: %d\n", n_teams))
cat(sprintf("Matches used: %d\n", nrow(df_pois)))

#Results & save
print(head(pois_abilities, 20))

top20_pois_table <- pois_abilities %>%
  slice_head(n = 20) %>%
  gt() %>%
  tab_header(
    title = "Top 20 Teams by Bivariate Poisson Net Strength",
    subtitle = "Goal-based ratings scaled to z-scores (2018-present)"
  ) %>%
  cols_label(
    rank_pois = "Rank",
    team = "Team",
    attack_strength = "Attack",
    defence_strength = "Defence",
    net_strength = "Net Strength"
  ) %>%
  fmt_number(
    columns = c(attack_strength, defence_strength, net_strength),
    decimals = 2
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = rank_pois, rows = rank_pois <= 3)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightblue"),
    locations = cells_body(rows = rank_pois <= 3)
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 12,
    table.width = pct(80)
  )

gtsave(top20_pois_table, "C:\\Users\\freya\\Downloads\\SMGT430\\top20_pois_table.pdf")

# =============================================================================
# SECTION 6: ELO RATING SYSTEM
# =============================================================================
# Features:
#   Competition-specific K-factor (World Cup = 60, Friendlies = 20)
#   Margin-of-victory multiplier: log(|gd| + 1) + 1
#   Draws handled as 0.5 score for both sides

elo_k_map <- c(
  "FIFA World Cup"                                 = 60,
  "CONMEBOL Copa America"                          = 50,
  "UEFA European Football Championship"            = 50,
  "Africa Cup of Nations"                          = 40,
  "AFC Asian Cup"                                  = 40,
  "CONCACAF Gold Cup"                              = 40,
  "FIFA World Cup Qualification - UEFA"            = 40,
  "FIFA World Cup Qualification - CONMEBOL"        = 40,
  "FIFA World Cup Qualification - AFC"             = 32,
  "FIFA World Cup Qualification - CAF"             = 32,
  "FIFA World Cup Qualification - CONCACAF"        = 32,
  "UEFA Nations League"                            = 35,
  "UEFA European Football Championship Qualifying" = 30,
  "International Friendlies (M)"                   = 20
)

elo_ratings <- setNames(rep(1500, length(all_teams)), all_teams)
df_elo_sorted <- df_clean %>%
  arrange(Date) %>%
  mutate(K = ifelse(Competition %in% names(elo_k_map),
                    elo_k_map[Competition], 30))

for (i in seq_len(nrow(df_elo_sorted))) {
  h <- df_elo_sorted$team1[i]; a <- df_elo_sorted$team2[i]
  if (!h %in% names(elo_ratings)) elo_ratings[h] <- 1500
  if (!a %in% names(elo_ratings)) elo_ratings[a] <- 1500
  Rh <- elo_ratings[h]; Ra <- elo_ratings[a]
  Eh <- 1 / (1 + 10^((Ra - Rh) / 400))
  g1 <- df_elo_sorted$score1[i]; g2 <- df_elo_sorted$score2[i]
  Sh <- case_when(g1 > g2 ~ 1, g1 < g2 ~ 0, TRUE ~ 0.5)
  mov_mult       <- log(abs(g1 - g2) + 1) + 1
  K              <- df_elo_sorted$K[i]
  elo_ratings[h] <- Rh + K * mov_mult * (Sh       - Eh)
  elo_ratings[a] <- Ra + K * mov_mult * ((1 - Sh) - (1 - Eh))
}

elo_df <- data.frame(
  team = names(elo_ratings),
  elo  = as.numeric(elo_ratings),
  stringsAsFactors = FALSE
) %>%
  filter(team %in% eligible_teams) %>%
  arrange(desc(elo)) %>%
  mutate(rank_elo = row_number())

#Results & save
print(head(elo_df, 20))

top20_elo_table <- elo_df %>%
  slice_head(n = 20) %>%
  gt() %>%
  tab_header(
    title = "Top 20 Teams by Elo Rating",
    subtitle = "Competition-weighted K-factor with margin-of-victory multiplier (2018-present)"
  ) %>%
  cols_label(
    rank_elo = "Rank",
    team = "Team",
    elo = "Elo Rating"
  ) %>%
  fmt_number(
    columns = c(elo),
    decimals = 0
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = rank_elo, rows = rank_elo <= 3)
  ) %>%
  tab_style(
    style = cell_fill(color = "lightgreen"),
    locations = cells_body(rows = rank_elo <= 3)
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 12,
    table.width = pct(80)
  )

gtsave(top20_elo_table, "C:\\Users\\freya\\Downloads\\SMGT430\\top20_elo_table.pdf")

# =============================================================================
# SECTION 7: COMBINED TOP 20 RANKINGS
# =============================================================================

# Find teams present in all three models
common_teams <- Reduce(intersect, list(
  bt_abilities$team,
  pois_abilities$team,
  elo_df$team
))
cat(sprintf("Teams present in all three models: %d\n", length(common_teams)))

# Merge all three rankings
combined <- bt_abilities %>%
  filter(team %in% common_teams) %>%
  dplyr::select(team, ability, rank_bt) %>%
  inner_join(
    pois_abilities %>%
      filter(team %in% common_teams) %>%
      dplyr::select(team, attack_strength, defence_strength, net_strength, rank_pois),
    by = "team"
  ) %>%
  inner_join(
    elo_df %>%
      filter(team %in% common_teams) %>%
      dplyr::select(team, elo, rank_elo),
    by = "team"
  )

# Calculate normalized scores and composite ranking
combined <- combined %>%
  mutate(
    norm_bt   = (ability - min(ability)) / (max(ability) - min(ability)),
    norm_pois = (net_strength - min(net_strength)) / (max(net_strength) - min(net_strength)),
    norm_elo  = (elo - min(elo)) / (max(elo) - min(elo)),
    composite = (norm_bt + norm_pois + norm_elo) / 3
  ) %>%
  arrange(desc(composite)) %>%
  mutate(rank_composite = row_number())

# Format combined data with 2 decimal places
combined_formatted <- combined %>%
  mutate(
    ability = round(ability, 2),
    net_strength = round(net_strength, 2),
    elo = round(elo, 0),  # Elo as whole number
    composite = round(composite, 4),
    norm_bt = round(norm_bt, 3),
    norm_pois = round(norm_pois, 3),
    norm_elo = round(norm_elo, 3)
  )

#Results & save
print(
  combined_formatted[1:20, c("rank_composite", "team", "ability", "net_strength", "elo", "composite")]
)

top20_overall <- combined_formatted %>%
  arrange(rank_composite) %>%
  slice_head(n = 20) %>%
  dplyr::select(
    rank_composite,
    team,
    ability,
    net_strength,
    elo,
    composite
  ) %>%
  gt() %>%
  tab_header(
    title = md("**Top 20 Overall Teams - 2026 FIFA World Cup**"),
    subtitle = md("*Composite ranking based on Bradley-Terry, Bivariate Poisson & Elo models*")
  ) %>%
  cols_label(
    rank_composite = md("**Rank**"),
    team = md("**Team**"),
    ability = md("**BT Ability**"),
    net_strength = md("**Poisson Net**"),
    elo = md("**Elo Rating**"),
    composite = md("**Composite Score**")
  ) %>%
  fmt_number(
    columns = c(ability, net_strength),
    decimals = 2
  ) %>%
  fmt_number(
    columns = c(elo),
    decimals = 0
  ) %>%
  fmt_number(
    columns = c(composite),
    decimals = 4
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold", color = "white"),
      cell_fill(color = "#1f77b4")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = "large"),
    locations = cells_body(columns = rank_composite, rows = rank_composite <= 3)
  ) %>%
  tab_style(
    style = cell_fill(color = "#ffd700"),  # Gold for top 3
    locations = cells_body(rows = rank_composite <= 3)
  ) %>%
  tab_style(
    style = cell_fill(color = "#c0c0c0"),  # Silver for ranks 4-5
    locations = cells_body(rows = rank_composite %in% 4:5)
  ) %>%
  tab_style(
    style = cell_fill(color = "#cd7f32"),  # Bronze for ranks 6-10
    locations = cells_body(rows = rank_composite %in% 6:10)
  ) %>%
  tab_footnote(
    footnote = "BT Ability: Time-weighted Bradley-Terry strength estimate (log scale)",
    locations = cells_column_labels(columns = ability)
  ) %>%
  tab_footnote(
    footnote = "Poisson Net: Attack minus defence strength from goal-based rating",
    locations = cells_column_labels(columns = net_strength)
  ) %>%
  tab_footnote(
    footnote = "Elo Rating: Competition-weighted with margin-of-victory multiplier",
    locations = cells_column_labels(columns = elo)
  ) %>%
  tab_footnote(
    footnote = "Composite Score: Normalized average of all three models (0-1 scale)",
    locations = cells_column_labels(columns = composite)
  ) %>%
  tab_source_note(
    source_note = "Data: 2018-present | Method: Equal weighting of BT, Poisson & Elo"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 18,
    heading.subtitle.font.size = 14,
    heading.title.font.weight = "bold",
    table.width = pct(90),
    table.background.color = "white",
    column_labels.background.color = "#f0f0f0",
    table.border.top.color = "#1f77b4",
    table.border.bottom.color = "#1f77b4",
    table.border.left.color = "white",
    table.border.right.color = "white"
  )

gtsave(top20_overall, "C:\\Users\\freya\\Downloads\\SMGT430\\top20_overall.pdf")

# =============================================================================
# SECTION 8: MONTE CARLO WORLD CUP SIMULATION (2026 FORMAT — 48 TEAMS)
# =============================================================================
# =============================================================================
# SECTION 8: MONTE CARLO WORLD CUP SIMULATION (2026 FORMAT — 48 TEAMS)
# =============================================================================
# Format: 48 teams, 12 groups of 4
# Group stage: top 2 from each group (24) + 8 best 3rd-place = 32 qualify
# Knockout: Round of 32 -> R16 -> QF -> SF -> Final
# =============================================================================

# Select top 48 teams for World Cup simulation
wc_teams_48 <- combined %>%
  arrange(desc(composite)) %>%
  slice_head(n = 48) %>%
  pull(team)

# Bradley-Terry win probability
bt_win_prob <- function(ti, tj) {
  if (!ti %in% bt_abilities$team || !tj %in% bt_abilities$team) return(0.5)
  li <- bt_abilities$ability[bt_abilities$team == ti]
  lj <- bt_abilities$ability[bt_abilities$team == tj]
  return(plogis(li - lj))
}

# Elo win probability (with draw probability)
elo_win_prob <- function(ti, tj) {
  Ri <- ifelse(ti %in% names(elo_ratings), elo_ratings[ti], 1500)
  Rj <- ifelse(tj %in% names(elo_ratings), elo_ratings[tj], 1500)
  pw <- 1 / (1 + 10^((Rj - Ri) / 400))
  pd <- 0.22 * 4 * pw * (1 - pw)  # Draw probability
  return(list(win = pw - pd/2, draw = pd, loss = 1 - pw - pd/2))
}

# Poisson match simulation (simplified for speed)
sim_pois_match <- function(team_h, team_a) {
  # Get net strengths
  ns_h <- pois_abilities$net_strength[pois_abilities$team == team_h]
  ns_a <- pois_abilities$net_strength[pois_abilities$team == team_a]
  if (length(ns_h) == 0) ns_h <- 0
  if (length(ns_a) == 0) ns_a <- 0
  # Expected goals (home advantage = 0.3 on log scale)
  lambda_h <- exp(0.3 + ns_h - ns_a)
  lambda_a <- exp(0.0 + ns_a - ns_h)
  # Small covariance for draws
  lambda_cov <- 0.15
  # Generate correlated Poisson goals
  z3 <- rpois(1, lambda_cov)
  return(list(g1 = rpois(1, lambda_h) + z3, g2 = rpois(1, lambda_a) + z3))
}

# Simulate a single match
simulate_match <- function(team_h, team_a, model, home_advantage = TRUE) {
  if (model == "bt") {
    p_win <- bt_win_prob(team_h, team_a)
    if (home_advantage) p_win <- p_win * 1.15 / (p_win * 1.15 + (1-p_win))
    r <- runif(1)
    if (r < p_win) return(list(g1 = 2, g2 = 1))  # Home win
    if (r < p_win + 0.22) return(list(g1 = 1, g2 = 1))  # Draw
    return(list(g1 = 1, g2 = 2))  # Away win
  }
  if (model == "pois") {
    res <- sim_pois_match(team_h, team_a)
    return(list(g1 = res$g1, g2 = res$g2))
  }
  if (model == "elo") {
    probs <- elo_win_prob(team_h, team_a)
    r <- runif(1)
    if (r < probs$win) return(list(g1 = 2, g2 = 1))
    if (r < probs$win + probs$draw) return(list(g1 = 1, g2 = 1))
    return(list(g1 = 1, g2 = 2))
  }
}

# Group Stage Simulation
simulate_group_stage <- function(group_teams, model) {
  n_teams <- length(group_teams)
  pts <- rep(0, n_teams)
  gd <- rep(0, n_teams)
  names(pts) <- group_teams
  names(gd) <- group_teams
  # Play each pair once
  for (i in 1:(n_teams-1)) {
    for (j in (i+1):n_teams) {
      team_h <- group_teams[i]
      team_a <- group_teams[j]
      # Home/away assigned randomly for fairness
      if (runif(1) < 0.5) {
        res <- simulate_match(team_h, team_a, model, home_advantage = TRUE)
        g1 <- res$g1; g2 <- res$g2
      } else {
        res <- simulate_match(team_a, team_h, model, home_advantage = TRUE)
        g1 <- res$g2; g2 <- res$g1
      }
      # Update points and goal difference
      if (g1 > g2) {
        pts[team_h] <- pts[team_h] + 3
      } else if (g2 > g1) {
        pts[team_a] <- pts[team_a] + 3
      } else {
        pts[team_h] <- pts[team_h] + 1
        pts[team_a] <- pts[team_a] + 1
      }
      gd[team_h] <- gd[team_h] + (g1 - g2)
      gd[team_a] <- gd[team_a] + (g2 - g1)
    }
  }
  # Return sorted standings
  standings <- data.frame(team = group_teams, pts = pts, gd = gd)
  standings <- standings[order(-standings$pts, -standings$gd), ]
  return(standings)
}

# Knockout Match Simulation
simulate_knockout <- function(team_h, team_a, model) {
  # Neutral venue, no home advantage
  res <- simulate_match(team_h, team_a, model, home_advantage = FALSE)
  if (res$g1 > res$g2) return(team_h)
  if (res$g2 > res$g1) return(team_a)
  # Penalty shootout (50/50 if draw after 90 mins)
  return(if (runif(1) < 0.5) team_h else team_a)
}

# Full Tournament Simulation
simulate_tournament <- function(model) {
  # Create 12 groups of 4 (seeded by composite rank)
  groups <- vector("list", 12)
  for (g in 1:12) {
    groups[[g]] <- wc_teams_48[seq(g, 48, by = 12)]
  }
  # Step 2: Group stage
  group_winners <- character()
  group_runners <- character()
  third_placed <- data.frame()
  for (g in 1:12) {
    standings <- simulate_group_stage(groups[[g]], model)
    group_winners <- c(group_winners, standings$team[1])
    group_runners <- c(group_runners, standings$team[2])
    
    third_placed <- rbind(third_placed, 
                          data.frame(team = standings$team[3], 
                                     pts = standings$pts[3], 
                                     gd = standings$gd[3]))
  }
  # Step 3: Select 8 best third-placed teams
  third_placed <- third_placed[order(-third_placed$pts, -third_placed$gd), ]
  best_thirds <- third_placed$team[1:8]
  # Step 4: Round of 32 (qualifiers)
  r32_teams <- c(group_winners, group_runners, best_thirds)
  r32_teams <- sample(r32_teams)  # Random draw
  # Step 5: Knockout rounds
  knockout_round <- function(teams) {
    winners <- character()
    for (i in seq(1, length(teams), by = 2)) {
      if (i + 1 <= length(teams)) {
        winner <- simulate_knockout(teams[i], teams[i+1], model)
        winners <- c(winners, winner)
      }
    }
    return(winners)
  }
  # R32 -> R16 -> QF -> SF -> Final
  r16 <- knockout_round(r32_teams)
  qf <- knockout_round(r16)
  sf <- knockout_round(qf)
  final_winner <- knockout_round(sf)
  
  return(final_winner)
}

# Run Simulations
n_simulations <- 10000
set.seed(2026)

# Run simulations for each model
run_model_sims <- function(model_name, model_code) {
  cat(sprintf("  Simulating %s model...", model_name))
  pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)
  champions <- character(n_simulations)
  for (i in 1:n_simulations) {
    champions[i] <- simulate_tournament(model_code)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # Calculate probabilities
  champ_table <- table(champions)
  champ_df <- data.frame(
    champions = names(champ_table),
    freq = as.numeric(champ_table),
    prob = as.numeric(champ_table) / n_simulations,
    model = model_name,
    stringsAsFactors = FALSE
  )
  champ_df <- champ_df[order(-champ_df$prob), ]
  return(champ_df)
}

# Run all three models
champ_bt <- run_model_sims("Bradley-Terry", "bt")
champ_pois <- run_model_sims("Bivariate Poisson", "pois")
champ_elo <- run_model_sims("Elo", "elo")

# Results
# Display top 10 for each model
for (df in list(champ_bt, champ_pois, champ_elo)) {
  model_name <- unique(df$model)
  cat(sprintf("【 %s MODEL 】\n", model_name))
  cat(strrep("-", 40), "\n")
  top10 <- head(df, 10)
  for (i in 1:nrow(top10)) {
    cat(sprintf("  %2d. %-20s %5.1f%%\n", 
                i, top10$champions[i], top10$prob[i] * 100))
  }
  cat("\n")
}

# Get top 10 from each model
bt_top10 <- champ_bt %>%
  slice_head(n = 10) %>%
  dplyr::select(Champion = champions, BT_Prob = prob)

pois_top10 <- champ_pois %>%
  slice_head(n = 10) %>%
  dplyr::select(Champion = champions, Poisson_Prob = prob)

elo_top10 <- champ_elo %>%
  slice_head(n = 10) %>%
  dplyr::select(Champion = champions, Elo_Prob = prob)

# Combine all top 10 teams (unique)
all_top_teams <- unique(c(bt_top10$Champion, pois_top10$Champion, elo_top10$Champion))

# Create full comparison table
comparison_table <- data.frame(Champion = all_top_teams) %>%
  left_join(bt_top10, by = "Champion") %>%
  left_join(pois_top10, by = "Champion") %>%
  left_join(elo_top10, by = "Champion") %>%
  mutate(
    BT_Prob = round(BT_Prob * 100, 1),
    Poisson_Prob = round(Poisson_Prob * 100, 1),
    Elo_Prob = round(Elo_Prob * 100, 1)
  ) %>%
  replace(is.na(.), 0) %>%
  arrange(desc(BT_Prob + Poisson_Prob + Elo_Prob)) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank, Champion, BT_Prob, Poisson_Prob, Elo_Prob)

# Add average probability column
comparison_table$Avg_Prob <- round(
  (comparison_table$BT_Prob + comparison_table$Poisson_Prob + comparison_table$Elo_Prob) / 3, 1
)

# Create styled gt table
montecarlo_table <- comparison_table %>%
  gt() %>%
  tab_header(
    title = md("**2026 FIFA World Cup - Monte Carlo Simulation Results**"),
    subtitle = md(paste0("*Top contenders across three prediction models | ", 
                         format(n_simulations, big.mark = ","), 
                         " simulations per model*"))
  ) %>%
  cols_label(
    Rank = md("**Rank**"),
    Champion = md("**Team**"),
    BT_Prob = md("**Bradley-Terry**<br><span style='font-size:10px'>(%)</span>"),
    Poisson_Prob = md("**Bivariate Poisson**<br><span style='font-size:10px'>(%)</span>"),
    Elo_Prob = md("**Elo**<br><span style='font-size:10px'>(%)</span>"),
    Avg_Prob = md("**Average**<br><span style='font-size:10px'>(%)</span>")
  ) %>%
  fmt_number(
    columns = c(BT_Prob, Poisson_Prob, Elo_Prob, Avg_Prob),
    decimals = 1,
    suffix = "%"
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold", color = "white"),
      cell_fill(color = "#1f77b4")
    ),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", size = "large"),
    locations = cells_body(columns = Rank, rows = Rank <= 3)
  ) %>%
  tab_style(
    style = cell_fill(color = "#ffd700"),  # Gold for top 3
    locations = cells_body(rows = Rank <= 3)
  ) %>%
  tab_style(
    style = cell_fill(color = "#c0c0c0"),  # Silver for ranks 4-5
    locations = cells_body(rows = Rank %in% 4:5)
  ) %>%
  tab_style(
    style = cell_fill(color = "#cd7f32"),  # Bronze for ranks 6-10
    locations = cells_body(rows = Rank %in% 6:10)
  ) %>%
  # Color scale for probabilities
  tab_style(
    style = cell_text(color = "green", weight = "bold"),
    locations = cells_body(columns = BT_Prob, rows = BT_Prob > 15)
  ) %>%
  tab_style(
    style = cell_text(color = "green", weight = "bold"),
    locations = cells_body(columns = Poisson_Prob, rows = Poisson_Prob > 15)
  ) %>%
  tab_style(
    style = cell_text(color = "green", weight = "bold"),
    locations = cells_body(columns = Elo_Prob, rows = Elo_Prob > 15)
  ) %>%
  tab_style(
    style = cell_text(color = "orange", weight = "bold"),
    locations = cells_body(columns = BT_Prob, rows = BT_Prob > 5 & BT_Prob <= 15)
  ) %>%
  tab_style(
    style = cell_text(color = "orange", weight = "bold"),
    locations = cells_body(columns = Poisson_Prob, rows = Poisson_Prob > 5 & Poisson_Prob <= 15)
  ) %>%
  tab_style(
    style = cell_text(color = "orange", weight = "bold"),
    locations = cells_body(columns = Elo_Prob, rows = Elo_Prob > 5 & Elo_Prob <= 15)
  ) %>%
  # Add summary statistics footer
  tab_footnote(
    footnote = md("**Bradley-Terry**: Time-weighted strength estimates with competition importance"),
    locations = cells_column_labels(columns = BT_Prob)
  ) %>%
  tab_footnote(
    footnote = md("**Bivariate Poisson**: Goal-based ratings with attack/defence decomposition"),
    locations = cells_column_labels(columns = Poisson_Prob)
  ) %>%
  tab_footnote(
    footnote = md("**Elo**: Competition-weighted K-factor with margin-of-victory multiplier"),
    locations = cells_column_labels(columns = Elo_Prob)
  ) %>%
  tab_source_note(
    source_note = md(paste0(
      "*Note*: Values represent percentage probability of winning the 2026 FIFA World Cup | ",
      "Based on 48-team format (12 groups of 4) | ",
      "Group stage: top 2 + 8 best 3rd-place advance to Round of 32 | ",
      "Knockout: Single-elimination with penalty shootouts for draws"
    ))
  ) %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 18,
    heading.subtitle.font.size = 13,
    heading.title.font.weight = "bold",
    table.width = pct(95),
    table.background.color = "white",
    column_labels.background.color = "#f0f0f0",
    table.border.top.color = "#1f77b4",
    table.border.bottom.color = "#1f77b4",
    table.border.left.color = "white",
    table.border.right.color = "white",
    data_row.padding = px(8)
  )

gtsave(montecarlo_table, "C:\\Users\\freya\\Downloads\\SMGT430\\montecarlo.pdf")



# =============================================================================
# SECTION 9: Visualization- Elo Rating Progression for Top 10 Teams
# =============================================================================
wc_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(size = 11, colour = "grey40"),
    plot.caption     = element_text(size = 9,  colour = "grey55"),
    axis.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

top10_elo_teams <- elo_df$team[1:10]
elo_track2      <- setNames(rep(1500, length(all_teams)), all_teams)
df_hist <- df_clean %>%
  arrange(Date) %>%
  mutate(
    K  = ifelse(Competition %in% names(elo_k_map), elo_k_map[Competition], 30),
    ym = floor_date(Date, "month")
  )
months_seq       <- sort(unique(df_hist$ym))
snapshot_ratings <- data.frame()

for (ym_int in as.integer(months_seq)) {
  ym_date <- as.Date(ym_int, origin = "1970-01-01")
  block   <- df_hist %>% filter(ym == ym_date)
  
  for (i in seq_len(nrow(block))) {
    h <- block$team1[i]; a <- block$team2[i]
    if (!h %in% names(elo_track2)) elo_track2[h] <- 1500
    if (!a %in% names(elo_track2)) elo_track2[a] <- 1500
    
    Rh <- elo_track2[h]; Ra <- elo_track2[a]
    Eh <- 1 / (1 + 10^((Ra - Rh) / 400))
    g1 <- block$score1[i]; g2 <- block$score2[i]
    Sh      <- case_when(g1 > g2 ~ 1, g1 < g2 ~ 0, TRUE ~ 0.5)
    mov_m   <- log(abs(g1 - g2) + 1) + 1
    K       <- block$K[i]
    elo_track2[h] <- Rh + K * mov_m * (Sh       - Eh)
    elo_track2[a] <- Ra + K * mov_m * ((1 - Sh) - (1 - Eh))
  }
  snap <- data.frame(
    date = ym_date,
    team = top10_elo_teams,
    elo  = elo_track2[top10_elo_teams]
  )
  snapshot_ratings <- rbind(snapshot_ratings, snap)
}

end_labels <- snapshot_ratings %>%
  group_by(team) %>%
  slice_tail(n = 1)

p4 <- ggplot(snapshot_ratings,
             aes(x = date, y = elo, colour = team, group = team)) +
  geom_line(linewidth = 1.05, alpha = 0.9) +
  geom_point(data = end_labels, size = 2.8) +
  geom_label_repel(data = end_labels, aes(label = team),
                   size = 3.1, nudge_x = 50,
                   direction = "y", segment.colour = "grey55") +
  scale_colour_viridis_d(option = "turbo", guide = "none") +
  scale_x_date(date_labels = "%b '%y", date_breaks = "6 months") +
  labs(
    title    = "Elo Rating Progression — Top 10 Teams (2018 to Present)",
    subtitle = "Competition-weighted K-factor | Margin-of-victory multiplier applied",
    x = NULL, y = "Elo Rating",
    caption  = "All teams initialised at 1500 | Ratings updated match-by-match in chronological order"
  ) +
  wc_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(p4)
ggsave("vis4_elo_progression.png", p4, width = 13, height = 7, dpi = 150)

# =============================================================================
# SECTION 10: Visualization- Model Agreement Heatmap
# =============================================================================
# Shows how well the three models agree with each other
# Calculate correlation between model rankings
rank_cor_matrix <- combined %>%
  dplyr::select(rank_bt, rank_pois, rank_elo) %>%
  cor(method = "spearman")

# Create heatmap
cor_df <- as.data.frame(rank_cor_matrix) %>%
  rownames_to_column("Model1") %>%
  pivot_longer(cols = -Model1, names_to = "Model2", values_to = "Correlation") %>%
  mutate(
    Model1 = recode(Model1, rank_bt = "BT", rank_pois = "Poisson", rank_elo = "Elo"),
    Model2 = recode(Model2, rank_bt = "BT", rank_pois = "Poisson", rank_elo = "Elo")
  )

p6 <- ggplot(cor_df, aes(x = Model1, y = Model2, fill = Correlation, label = round(Correlation, 2))) +
  geom_tile(color = "white", size = 1) +
  geom_text(size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#E63946", mid = "white", high = "#2A9D8F", 
                       midpoint = 0.5, limits = c(0, 1),
                       name = "Spearman\nCorrelation") +
  labs(
    title = "Model Agreement: Ranking Correlation Heatmap",
    subtitle = "Higher correlation indicates stronger agreement between models",
    x = NULL, y = NULL,
    caption = "Based on Spearman rank correlation of top 48 teams"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(p6)
ggsave("vis6_model_agreement_heatmap.png", p6, width = 8, height = 6, dpi = 150)

# -----------------------------------------------------------------------------
# SECTION 11: Visualization- Champion Probability Distribution - Top Contenders Race
# -----------------------------------------------------------------------------
# Shows the relative probabilities and creates a "horseshoe" race plot
# Prepare champion probability data for top 8 teams across models
top8_teams <- unique(c(champ_bt$champions[1:8], 
                       champ_pois$champions[1:8], 
                       champ_elo$champions[1:8]))

champion_race <- all_champs %>%
  filter(champions %in% top8_teams) %>%
  mutate(
    champions = factor(champions, levels = rev(top8_teams)),
    model = factor(model, levels = c("Bradley-Terry", "Bivariate Poisson", "Elo"))
  )

p7 <- ggplot(champion_race, aes(x = prob, y = champions, fill = model)) +
  geom_col(position = position_dodge(0.85), width = 0.75, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", prob * 100)),
            position = position_dodge(0.85), hjust = -0.15, size = 3.2) +
  scale_fill_manual(values = c("Bradley-Terry" = "#E63946",
                               "Bivariate Poisson" = "#457B9D",
                               "Elo" = "#2A9D8F"),
                    name = "Prediction Model") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, max(champion_race$prob) * 1.15)) +
  labs(
    title = "World Cup Champion Probabilities: Top Contenders",
    subtitle = "Comparison across three prediction models (10,000 simulations each)",
    x = "Probability of Winning Tournament",
    y = NULL,
    caption = "Error bars show model uncertainty | Based on 2026 format (48 teams)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p7)
ggsave("vis7_champion_race.png", p7, width = 11, height = 7, dpi = 150)


# -----------------------------------------------------------------------------
# SECTION 12: Visualization- Model Uncertainty - Standard Error Comparison
# -----------------------------------------------------------------------------
# Shows which teams have the most uncertain ratings
# Prepare uncertainty data
uncertainty_data <- bt_abilities %>%
  dplyr::select(team, ability, se_bt) %>%
  mutate(
    ci_lower = ability - 1.96 * se_bt,
    ci_upper = ability + 1.96 * se_bt,
    uncertainty = se_bt
  ) %>%
  arrange(desc(uncertainty)) %>%
  slice_head(n = 20)

p9 <- ggplot(uncertainty_data, aes(x = reorder(team, uncertainty), y = uncertainty)) +
  geom_col(aes(fill = uncertainty), width = 0.7) +
  geom_text(aes(label = round(uncertainty, 2)), hjust = -0.2, size = 3.2) +
  scale_fill_gradient(low = "#2A9D8F", high = "#E63946", 
                      name = "Standard Error") +
  coord_flip() +
  labs(
    title = "Model Uncertainty: Teams with Widest Confidence Intervals",
    subtitle = "Higher standard error = less precise Bradley-Terry ability estimate",
    x = NULL,
    y = "Standard Error (BT Ability)",
    caption = "Based on time-weighted Bradley-Terry model | 95% CI = ability ± 1.96×SE"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

print(p9)
ggsave("vis9_model_uncertainty.png", p9, width = 10, height = 7, dpi = 150)
cat("  Saved vis9_model_uncertainty.png\n")