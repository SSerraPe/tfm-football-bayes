# Football Factor Interpretation — La Liga Low-Rank Model

Model: low-rank Sigma_a (Q_a = 2) + diagonal Sigma_b
Data: real La Liga Wyscout data (1650 players, 12 seasons, 53 features)

## Inferred Factor Labels

### Factor_1: Attacking / Creative

**Top positive loadings:** per90_forward_passes (-0.816), per90_recoveries (-0.737), per90_progressive_passes (-0.694), per90_touch_in_box (0.691), per90_offensive_duels (0.671)
**Top negative loadings:** per90_forward_passes (-0.816), per90_recoveries (-0.737), per90_progressive_passes (-0.694), per90_lateral_passes (-0.658), per90_passes (-0.629)

**Interpretation:**
Factor 1 captures the latent dimension labelled 'Attacking / Creative'.
Players with high scores on this factor tend to show elevated values
on the positively-loading features and suppressed values on negatively-loading ones.

### Factor_2: Attacking / Creative

**Top positive loadings:** per90_assists (0.297), per90_forward_passes (-0.214), per90_recoveries (-0.204), per90_offensive_duels (0.202), per90_touch_in_box (0.199)
**Top negative loadings:** per90_forward_passes (-0.214), per90_recoveries (-0.204), per90_progressive_passes (-0.176), per90_long_passes (-0.171), per90_interceptions (-0.168)

**Interpretation:**
Factor 2 captures the latent dimension labelled 'Attacking / Creative'.
Players with high scores on this factor tend to show elevated values
on the positively-loading features and suppressed values on negatively-loading ones.

## Football Context

In modern football analytics, player profiles can be decomposed into a small number
of latent axes that capture correlated patterns of action. The two-factor structure
found here is consistent with the literature distinguishing:

- **Attacking output** players: frequent goal contributions, touches in box, progressive runs,
  dribbles, shot-creating actions.
- **Defensive / physical** players: high recovery counts, pressing duels, defensive duels,
  interceptions, clearances.

The factor loadings have sign (+/-) which indicates direction:
- A positive loading means higher factor score => higher value on that feature.
- A negative loading means higher factor score => lower value on that feature.

## Caveats

1. Anchor constraint: features 1 and 2 are set as positive anchors for factors 1 and 2
   respectively. This fixes the sign convention for identification.
2. Season effects are modelled as diagonal (Sigma_b = Diag(sigma_b^2)). Any
   cross-feature seasonal co-movement is absorbed into player or residual components.
3. With Q_a = 2, the model captures the two dominant axes of player variation.
   A third factor may exist but requires more data for reliable estimation.

## Summary

Posterior summary: /Users/sserra/Documents/01_MESIO_UPC/TFM/model/outputs/tables/10_real_lowrank_a_diag_b_posterior_summary.csv
Loading table:     /Users/sserra/Documents/01_MESIO_UPC/TFM/model/outputs/tables/11_real_lambda_a_loading_table.csv
Player scores:     /Users/sserra/Documents/01_MESIO_UPC/TFM/model/outputs/tables/12_player_factor_scores.csv

Generated: 2026-05-29 01:40:21.517533
