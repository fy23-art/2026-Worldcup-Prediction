# 2026-Worldcup-Prediction

International football lacks the high-frequency match data found in domestic leagues, making "current form" difficult to quantify. Accurate forecasting allows national federations to optimize resource allocation, such as prioritizing friendlies that maximize Elo gains for better tournament seeding. Furthermore, sponsors and broadcasters use these insights for strategic planning based on the probabilistic likelihood of specific teams advancing to high-value knockout stages. 


To evaluate the competitive landscape of the 2026 FIFA World Cup, this study benchmarks a time-weighted Bradley-Terry model, integrated within a bivariate Poisson framework, against traditional Elo rating systems. By assessing the top 20 nations by latent ability, we aim to determine which modeling approach most effectively captures the winning probability. The following sections detail the primary methodologies employed to derive these probabilistic forecasts:


Primary data for this analysis was programmatically extracted using the worldfootballR library. To ensure the integrity and predictive relevance of the dataset, the following processing steps and feature engineering techniques were applied:
- Team Identification: Standardized mapping of team names across multiple seasons and competitions.
- Venue Assignment: Categorization of fixtures into Home and Away contexts to account for territorial advantage.
- Competition Weighting: Implementation of importance coefficients based on the tier of competition. Matches were weighted hierarchically (e.g., FIFA World Cup and UEFA European Championship fixtures carried higher statistical significance than International Friendlies).
- Temporal Decay: Application of a time-decay function to the historical data, ensuring that recent match outcomes have a greater influence on the model than older results.
- Sparsity Filtering: To maintain statistical significance and mitigate noise from outliers, a minimum threshold was established. Teams with fewer than 8 total matches were excluded from the final dataset.
