# Validation Interpretation

Generated: 2026-05-14

## Execution Status

All three fitted models are present and usable through lightweight fit pointers plus persistent CmdStan CSVs:

- `03_real_diagonal_additive`
- `05_sim_diagonal_additive`
- `06_sim_lowrank_recovery`

The original `07` memory failure was caused by optional trace-plot draw extraction from very large CmdStan CSVs. Stage `07` was patched to use compact posterior summaries and no longer loads full draw matrices. Stage `08` was also patched to skip LOO unless `BFA_RUN_LOO=true`.

## Sampler Geometry

All three models have clean HMC geometry:

| model | draws | divergences | max treedepth | treedepth hits | min E-BFMI | mean accept |
|---|---:|---:|---:|---:|---:|---:|
| real diagonal | 4000 | 0 | 8 | 0 | 0.561 | 0.953 |
| simulated diagonal | 4000 | 0 | 8 | 0 | 0.551 | 0.951 |
| simulated low-rank | 4000 | 0 | 8 | 0 | 0.586 | 0.951 |

There is no evidence of divergent transitions or treedepth saturation. The main issue is posterior mixing/identifiability, not numerical HMC failure.

## Posterior Health

| model | summarized parameters | max Rhat | Rhat > 1.01 | Rhat > 1.05 | min bulk ESS | min tail ESS |
|---|---:|---:|---:|---:|---:|---:|
| real diagonal | 531 | 1.035 | 13 | 0 | 159 | 175 |
| simulated diagonal | 531 | 1.010 | 3 | 0 | 779 | 1412 |
| simulated low-rank | 688 | 1.739 | 233 | 208 | 6 | 107 |

The real diagonal model is mostly usable but has mild convergence weakness in a few player-variance parameters. The simulated diagonal model is healthy. The simulated low-rank recovery model is not converged for direct loading interpretation.

## Real-Data Diagonal Model

The real diagonal model estimates most variation as player-level for volume/intensity features, with smaller season and residual components.

Highest player-level shares:

- `per90_recoveries`: 88.4% player, 1.6% season, 10.0% residual
- `per90_touch_in_box`: 87.2% player, 0.1% season, 12.7% residual
- `per90_offensive_duels`: 86.9% player, 5.6% season, 7.5% residual
- `per90_shots`: 86.0% player, 0.8% season, 13.2% residual
- `per90_forward_passes`: 84.6% player, 0.4% season, 15.0% residual

Highest season-level shares:

- `per90_smart_passes`: 44.0% player, 37.6% season, 18.3% residual
- `per90_pressing_duels`: 45.5% player, 32.0% season, 22.5% residual
- `per90_accelerations`: 59.8% player, 16.4% season, 23.8% residual
- `per90_dribbles_against`: 69.3% player, 13.6% season, 17.1% residual

Highest residual shares:

- `rate_successful_crosses`: 97.5% residual
- `rate_successful_sliding_tackles`: 93.2% residual
- `rate_shots_on_target`: 83.2% residual
- `rate_goal_conversion`: 81.8% residual
- `rate_dribbles_against_won`: 80.8% residual

Interpretation: stable player differences dominate many event-volume features. Percent/rate variables are noisier and much less attributable to stable player or season effects.

## Simulated Diagonal Fit To Low-Rank Data

The diagonal model fitted to simulated low-rank data recovers covariance diagonals well overall:

| component | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|
| player | 0.0078 | 0.0112 | 0.996 | 0.906 |
| residual | 0.0040 | 0.0055 | 0.998 | 0.925 |
| season | 0.0300 | 0.0468 | 0.885 | 0.868 |

Player and residual diagonal recovery are excellent. Season variance recovery is weaker and biased upward, which is expected because there are only 12 seasons and the low-rank season signal is harder to identify from the observed panel.

## Simulated Low-Rank Recovery

The low-rank recovery model has clean sampler geometry but poor posterior mixing for loading parameters:

| parameter group | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|
| `Lambda_a` | 0.082 | 0.131 | 0.804 | 0.840 |
| `Lambda_b` | 0.102 | 0.124 | 0.262 | 1.000 |
| `Sigma_a_diag` | 0.009 | 0.013 | 0.997 | 0.811 |
| `Sigma_b_diag` | 0.057 | 0.093 | 0.886 | 0.585 |
| `psi_a_sd` | 0.011 | 0.022 | 0.941 | 0.925 |
| `psi_b_sd` | 0.040 | 0.052 | 0.740 | 0.925 |
| `sigma_epsilon_sd` | 0.004 | 0.005 | 0.998 | 0.925 |

The model recovers residual scales and player covariance diagonals well, but loading-level recovery is not reliable. The worst Rhat values are around 1.74 for `lambda_a_diag[1]`, `lambda_b_diag[1]`, and many free loadings. This is consistent with weak identification/multimodality between factor loadings and uniqueness terms, especially for the season block.

Important detail: the low-rank model is not failing because of divergences. It is failing because the posterior has poorly mixed loading/uniqueness modes under the current constraints and data size.

## PSIS-LOO

PSIS-LOO was computed after installing the missing `loo` package. The output is saved in:

- `outputs/tables/08_loo_summary.csv`

| model | elpd_loo | se_elpd_loo | p_loo | looic | Pareto k > 0.7 | Pareto k > 1 | max Pareto k |
|---|---:|---:|---:|---:|---:|---:|---:|
| real diagonal | -224673.4 | 1007.0 | 49545.7 | 449346.9 | 4354 | 3088 | 3.82 |
| simulated diagonal | -201763.1 | 394.7 | 41871.8 | 403526.1 | 4529 | 2878 | 2.88 |
| simulated low-rank | -196186.8 | 363.3 | 32963.7 | 392373.7 | 4215 | 2104 | 2.08 |

The simulated low-rank model has higher expected log predictive density than the simulated diagonal model on the simulated data, as expected because the data were generated from a low-rank process. The nominal difference is large:

\[
\Delta \mathrm{elpd}_{\mathrm{lowrank-diagonal}} \approx 5576.
\]

However, these PSIS-LOO estimates are not reliable as formal model-selection evidence because almost all observations have problematic Pareto-k values. In particular, all three models have thousands of observations with \(k > 0.7\), and thousands with \(k > 1\). This means the importance-sampling approximation is unstable. The LOO numbers are useful as a warning about influential observations and possible overfitting/high leverage, but not as final model-comparison evidence.

For defensible predictive comparison, use exact K-fold cross-validation or refit-based leave-one-out for a targeted subset. Moment matching may help, but with this many high Pareto-k observations, K-fold is the cleaner next step.

## Acceptance Decision

- Accept the real diagonal model as a defensible descriptive additive variance decomposition, with caveats for a small number of Rhat/ESS warnings.
- Accept the simulated diagonal model as a successful sanity check for diagonal variance recovery under misspecified low-rank data.
- Do not accept the current low-rank recovery model for direct loading interpretation.
- Treat low-rank covariance diagonal recovery as partially informative for player effects, but not sufficient for reliable loading recovery or season-block interpretation.
- Do not use the current PSIS-LOO estimates as decisive model-selection evidence because Pareto-k diagnostics fail badly.

## Recommended Next Steps

1. Keep using the diagonal additive model for the thesis-level real-data decomposition unless the research question requires loading-level recovery.
2. If low-rank recovery is required, refit `06` with stronger identifiability:
   - initialize near the simulated truth for recovery experiments,
   - use stronger positive anchor priors,
   - consider lower-bounding anchor loadings away from zero,
   - add stronger shrinkage or ordering constraints,
   - consider marginalizing latent effects if computationally feasible.
3. Treat season low-rank recovery as inherently fragile with only 12 seasons.
4. Avoid `BFA_RUN_LOO=true` until there is a specific need, because log-likelihood draw extraction can again be memory-heavy.
