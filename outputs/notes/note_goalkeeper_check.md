# Goalkeeper Structural Check

## Finding

Goalkeepers were **not** excluded from the model panel. Only GK-specific
feature columns (prefix `gk_` or `goal_kicks`) were filtered in stage 01.
GK flag: `total_gk_saves > 0` joined from the raw Wyscout CSV by player_id/season_id.

- **358 GK player-season rows** out of N=4944 (7.2%)
- **121 distinct GK players** present in `model_objects$metadata`

## GK-confounded features (|icc_all - icc_outfield| > 0.02)

- `rate_defensive_duels_won`: icc_all=0.243, icc_outfield=0.389 (delta=-0.146)
- `per90_own_half_losses`: icc_all=0.671, icc_outfield=0.528 (delta=0.143)
- `per90_defensive_duels`: icc_all=0.815, icc_outfield=0.709 (delta=0.106)
- `per90_loose_ball_duels`: icc_all=0.771, icc_outfield=0.677 (delta=0.095)
- `rate_aerial_duels_won`: icc_all=0.768, icc_outfield=0.673 (delta=0.095)
- `per90_dangerous_opponent_half_recoveries`: icc_all=0.475, icc_outfield=0.382 (delta=0.093)
- `per90_fouls`: icc_all=0.732, icc_outfield=0.645 (delta=0.087)
- `per90_losses`: icc_all=0.785, icc_outfield=0.702 (delta=0.084)
- `per90_opponent_half_recoveries`: icc_all=0.792, icc_outfield=0.710 (delta=0.083)
- `rate_dribbles_against_won`: icc_all=0.219, icc_outfield=0.300 (delta=-0.081)
- `per90_dribbles_against`: icc_all=0.756, icc_outfield=0.695 (delta=0.061)
- `per90_back_passes`: icc_all=0.830, icc_outfield=0.770 (delta=0.059)
- `rate_offensive_duels_won`: icc_all=0.416, icc_outfield=0.470 (delta=-0.054)
- `per90_yellow_cards`: icc_all=0.446, icc_outfield=0.396 (delta=0.050)
- `per90_dangerous_own_half_losses`: icc_all=0.482, icc_outfield=0.440 (delta=0.043)
- `per90_fouls_suffered`: icc_all=0.782, icc_outfield=0.746 (delta=0.036)
- `per90_key_passes`: icc_all=0.660, icc_outfield=0.630 (delta=0.029)
- `per90_pressing_duels`: icc_all=0.553, icc_outfield=0.524 (delta=0.029)
- `per90_xg_assist`: icc_all=0.675, icc_outfield=0.648 (delta=0.027)
- `per90_assists`: icc_all=0.398, icc_outfield=0.372 (delta=0.026)
- `per90_sliding_tackles`: icc_all=0.661, icc_outfield=0.635 (delta=0.026)
- `per90_accelerations`: icc_all=0.643, icc_outfield=0.618 (delta=0.025)
- `rate_successful_passes`: icc_all=0.769, icc_outfield=0.794 (delta=-0.025)
- `per90_shots_blocked`: icc_all=0.661, icc_outfield=0.636 (delta=0.025)
- `per90_shot_assists`: icc_all=0.771, icc_outfield=0.746 (delta=0.025)
- `per90_progressive_run`: icc_all=0.778, icc_outfield=0.754 (delta=0.024)
- `per90_interceptions`: icc_all=0.815, icc_outfield=0.793 (delta=0.021)
- `per90_through_passes`: icc_all=0.668, icc_outfield=0.647 (delta=0.021)
- `per90_aerial_duels`: icc_all=0.836, icc_outfield=0.816 (delta=0.021)

## Impact on borderline features (model-based icc_player 0.20–0.35)

- `rate_defensive_duels_won`: icc_outfield=0.389 → still borderline
- `rate_successful_dribbles`: icc_outfield=0.245 → still borderline

## Recommended fix

Add a row-level GK filter in `scripts/01_read_clean_feature_engineering.R`
before the `eligible` step, using the `total_gk_saves` column from the
longitudinal panel produced by stage 00 (`data/interim/longitudinal_player_seasons.csv`).

```r
eligible <- panel |>
  filter(minutes >= pipeline$minimum_minutes_model,
         total_gk_saves == 0 | is.na(total_gk_saves))  # exclude GK rows
```

This goes in stage 01's rebuild, bundled with Task K's feature drops and transforms.

_Written by scripts/25_gk_icc_analysis.R on 2026-06-30_
