# Football Factor Interpretation — La Liga Low-Rank Model

Model: low-rank Sigma_a (Q_a = 2) + diagonal Sigma_b
Data: real La Liga Wyscout data (1650 players, 12 seasons, 53 features)

## Inferred Factor Labels

### Factor_1: Attacking / Creative

**Top positive loadings:** per90_offensive_duels (0.847), per90_touch_in_box (0.804), per90_shots (0.779), per90_forward_passes (-0.760), per90_recoveries (-0.739)
**Top negative loadings:** per90_forward_passes (-0.760), per90_recoveries (-0.739), per90_long_passes (-0.697), per90_clearances (-0.641), per90_progressive_passes (-0.623)

**Interpretation:**
Factor 1 captures the latent dimension labelled 'Attacking / Creative'.
Players with high scores on this factor tend to show elevated values
on the positively-loading features and suppressed values on negatively-loading ones.

### Factor_2: Attacking / Creative

**Top positive loadings:** per90_received_pass (0.797), per90_passes (0.796), per90_back_passes (0.679), per90_lateral_passes (0.635), per90_passes_to_final_third (0.612)
**Top negative loadings:** rate_aerial_duels_won (-0.337), per90_head_shots (-0.243), per90_aerial_duels (-0.210), per90_long_passes (-0.061), per90_xg_shot (-0.045)

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

Generated: 2026-05-29 08:22:13.050647
