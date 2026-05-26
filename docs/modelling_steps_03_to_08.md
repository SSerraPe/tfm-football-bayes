# Modelling Steps 03 To 08: From Real-Data Fit To Validation

This note explains the modelling stages after the model objects were prepared. The aim of these stages is to move from a real-data Bayesian fit to simulation-based checks and final validation. The modelling unit is one observed player-season row, with 53 scaled football features, 4,944 rows, 1,650 players, and 12 seasons.

## Step 03 - Fit The Real-Data Diagonal Additive Model

**What was done**

The real scaled response matrix was fitted with the additive model:

```text
y_ij = mu + player_i + season_j + residual_ij
```

The player, season, and residual covariance matrices were assumed diagonal. This means each feature gets its own player variance, season variance, and residual variance, but the model does not estimate cross-feature covariance in this real-data baseline.

**Why it was done this way**

This is the cleanest first Bayesian decomposition of the data. It answers the main descriptive question directly: for each football feature, how much variation is stable player-level variation, how much is common season-level variation, and how much remains unexplained at player-season level? The diagonal assumption keeps the model identifiable and computationally manageable with many players and 53 response features.

**Results**

Sampler health was good:

| metric | result |
|---|---:|
| posterior draws | 4,000 |
| divergences | 0 |
| max treedepth | 8 |
| max-treedepth hits | 0 |
| min E-BFMI | 0.561 |
| mean acceptance statistic | 0.953 |

Posterior convergence was mostly acceptable, with mild warnings:

| metric | result |
|---|---:|
| summarized parameters | 531 |
| max Rhat | 1.035 |
| parameters with Rhat > 1.01 | 13 |
| parameters with Rhat > 1.05 | 0 |
| min bulk ESS | 159 |
| min tail ESS | 175 |

The substantive result is that many volume/intensity features are mostly player-level:

| feature | player share | season share | residual share |
|---|---:|---:|---:|
| `per90_recoveries` | 88.4% | 1.6% | 10.0% |
| `per90_touch_in_box` | 87.2% | 0.1% | 12.7% |
| `per90_offensive_duels` | 86.9% | 5.6% | 7.5% |
| `per90_shots` | 86.0% | 0.8% | 13.2% |
| `per90_forward_passes` | 84.6% | 0.4% | 15.0% |

The highest season-level shares were:

| feature | player share | season share | residual share |
|---|---:|---:|---:|
| `per90_smart_passes` | 44.0% | 37.6% | 18.3% |
| `per90_pressing_duels` | 45.5% | 32.0% | 22.5% |
| `per90_accelerations` | 59.8% | 16.4% | 23.8% |
| `per90_dribbles_against` | 69.3% | 13.6% | 17.1% |

Rate features were mostly residual/noisy:

| feature | residual share |
|---|---:|
| `rate_successful_crosses` | 97.5% |
| `rate_successful_sliding_tackles` | 93.2% |
| `rate_shots_on_target` | 83.2% |
| `rate_goal_conversion` | 81.8% |
| `rate_dribbles_against_won` | 80.8% |

Main output files:

- `fits/03_real_diagonal_additive_fit.rds`
- `outputs/tables/03_real_diagonal_additive_posterior_summary.csv`
- `outputs/diagnostics/03_real_diagonal_additive_diagnostic_summary.csv`
- `outputs/tables/07_real_variance_decomposition.csv`

## Step 04 - Simulate Low-Rank Additive Data

**What was done**

A synthetic dataset was generated on the same observed player-season structure as the real data. The simulated data used a low-rank plus diagonal additive design:

- player effect rank: 2
- season effect rank: 1
- same number of rows, features, players, and seasons as the real modelling data

The pipeline saved the simulated response matrix, the true loadings, true variance components, and Stan data for both the diagonal and low-rank models.

**Why it was done this way**

The real data do not have known truth. A simulation gives known truth while preserving the real panel structure, including the number of players, number of seasons, missing player-seasons, and feature dimension. This makes validation more realistic than a toy simulation, because the recovery problem has the same irregular panel shape as the real analysis.

**Results**

| item | value |
|---|---:|
| simulation seed | 20260913 |
| rows `N` | 4,944 |
| features `P` | 53 |
| players | 1,650 |
| seasons | 12 |
| player rank | 2 |
| season rank | 1 |

Main output files:

- `data/processed/simulated_Y_scaled.csv`
- `data/processed/simulated_lowrank_truth.rds`
- `outputs/tables/04_simulation_truth_loadings.csv`
- `outputs/tables/04_simulation_truth_variance_components.csv`
- `outputs/tables/04_simulation_summary.csv`

## Step 05 - Fit The Diagonal Model To Simulated Low-Rank Data

**What was done**

The same diagonal additive model from step 03 was fitted to the simulated low-rank data.

**Why it was done this way**

This is a misspecification check. Because the simulated data were generated with low-rank covariance, the diagonal model is intentionally too simple. The point is to test whether the simpler real-data model can still recover the main diagonal variance components, even when the true covariance contains cross-feature structure.

**Results**

Sampler health was good:

| metric | result |
|---|---:|
| posterior draws | 4,000 |
| divergences | 0 |
| max treedepth | 8 |
| max-treedepth hits | 0 |
| min E-BFMI | 0.551 |
| mean acceptance statistic | 0.951 |

Posterior convergence was strong:

| metric | result |
|---|---:|
| summarized parameters | 531 |
| max Rhat | 1.010 |
| parameters with Rhat > 1.01 | 3 |
| parameters with Rhat > 1.05 | 0 |
| min bulk ESS | 779 |
| min tail ESS | 1,412 |

Recovery of the known variance diagonals:

| component | bias | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|---:|
| player | 0.0028 | 0.0078 | 0.0112 | 0.996 | 0.906 |
| residual | -0.0009 | 0.0040 | 0.0055 | 0.998 | 0.925 |
| season | 0.0250 | 0.0300 | 0.0468 | 0.885 | 0.868 |

The diagonal model recovered player and residual variances very well. Season variance was weaker and biased upward, which is expected because there are only 12 seasons and season effects are harder to identify.

Main output files:

- `fits/05_sim_diagonal_additive_fit.rds`
- `outputs/tables/05_sim_diagonal_additive_posterior_summary.csv`
- `outputs/tables/08_sim_diagonal_recovery_metrics.csv`
- `outputs/tables/08_sim_diagonal_recovery_by_feature.csv`

## Step 06 - Fit The Low-Rank Recovery Model To Simulated Data

**What was done**

A low-rank recovery model was fitted to the same simulated data. Unlike step 05, this model matched the data-generating structure: low-rank plus diagonal player and season covariance.

**Why it was done this way**

This is the direct recovery test. If the low-rank model is correctly specified and identifiable, it should recover the true loadings, uniqueness terms, residual scales, and covariance diagonals. This stage checks whether a richer covariance model is actually reliable enough to interpret.

**Results**

Sampler geometry was clean:

| metric | result |
|---|---:|
| posterior draws | 4,000 |
| divergences | 0 |
| max treedepth | 8 |
| max-treedepth hits | 0 |
| min E-BFMI | 0.586 |
| mean acceptance statistic | 0.951 |

But posterior mixing was poor:

| metric | result |
|---|---:|
| summarized parameters | 688 |
| max Rhat | 1.739 |
| parameters with Rhat > 1.01 | 233 |
| parameters with Rhat > 1.05 | 208 |
| min bulk ESS | 6 |
| min tail ESS | 107 |

Loading recovery:

| parameter | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|
| `Lambda_a` | 0.082 | 0.131 | 0.804 | 0.840 |
| `Lambda_b` | 0.102 | 0.124 | 0.262 | 1.000 |

Variance and scale recovery:

| parameter | MAE | RMSE | correlation | 90% coverage |
|---|---:|---:|---:|---:|
| `Sigma_a_diag` | 0.009 | 0.013 | 0.997 | 0.811 |
| `Sigma_b_diag` | 0.057 | 0.093 | 0.886 | 0.585 |
| `psi_a_sd` | 0.011 | 0.022 | 0.941 | 0.925 |
| `psi_b_sd` | 0.040 | 0.052 | 0.740 | 0.925 |
| `sigma_epsilon_sd` | 0.004 | 0.005 | 0.998 | 0.925 |

The model recovered residual scales and player covariance diagonals well, but loading-level recovery was not reliable. The issue is not HMC failure, because there were no divergences. The issue is weak identification and poor posterior mixing for loadings and uniqueness terms, especially in the season block.

Main output files:

- `fits/06_sim_lowrank_recovery_fit.rds`
- `outputs/tables/06_sim_lowrank_recovery_posterior_summary.csv`
- `outputs/tables/08_sim_lowrank_loading_recovery_metrics.csv`
- `outputs/tables/08_sim_lowrank_variance_recovery_metrics.csv`

## Step 07 - Produce Diagnostics, Recovery Tables, And Figures

**What was done**

The saved fits and simulation truth were converted into compact diagnostic tables and figures. This stage wrote:

- real-data empirical and posterior variance decomposition tables
- simulated diagonal variance-vs-truth tables
- simulated low-rank loading and variance recovery tables
- PNG figures for visual checks

Trace plots were skipped intentionally to avoid loading very large draw matrices from the CmdStan CSV files.

**Why it was done this way**

The modelling results need to be inspectable without repeatedly loading full posterior draws into memory. Compact summaries are enough for the main thesis diagnostics: variance shares, convergence summaries, recovery against truth, and visual posterior-mean-vs-truth comparisons.

**Results**

All required summaries were available:

| artifact | available |
|---|---:|
| real diagonal summary | TRUE |
| simulated diagonal summary | TRUE |
| simulated low-rank summary | TRUE |

Figure status:

| plot group | status |
|---|---|
| diagnostic plots | summary plots written; trace plots skipped to avoid large draw loading |

Main output files:

- `outputs/tables/07_real_variance_decomposition.csv`
- `outputs/tables/07_real_empirical_variance_decomposition.csv`
- `outputs/tables/07_sim_diagonal_variance_vs_truth.csv`
- `outputs/tables/07_sim_lowrank_recovery_loadings.csv`
- `outputs/tables/07_sim_lowrank_recovery_variances.csv`
- `outputs/figures/07_real_empirical_variance_decomposition.png`
- `outputs/figures/07_sim_diagonal_covariance_diag_vs_truth.png`
- `outputs/figures/07_sim_lowrank_loadings_vs_truth.png`
- `outputs/figures/07_sim_lowrank_variances_vs_truth.png`

## Step 08 - Validate All Models

**What was done**

The final validation stage checked the availability of all fits and summaries, sampler health, posterior health, simulation recovery, and optional PSIS-LOO summaries.

**Why it was done this way**

Validation combines three different questions that should not be confused:

1. Did HMC sample the posterior geometry cleanly?
2. Did the posterior chains mix well enough for the parameters we want to interpret?
3. On simulated data, did posterior estimates recover known truth?

A model can pass one check and fail another. That happened here: the low-rank model had clean sampler geometry but poor loading convergence.

**Results**

Artifact status:

| item | result |
|---|---:|
| fit files present | 3 / 3 |
| posterior summaries present | 3 / 3 |

Posterior health:

| model | parameters | max Rhat | Rhat > 1.01 | Rhat > 1.05 | min bulk ESS | min tail ESS |
|---|---:|---:|---:|---:|---:|---:|
| real diagonal | 531 | 1.035 | 13 | 0 | 159 | 175 |
| simulated diagonal | 531 | 1.010 | 3 | 0 | 779 | 1,412 |
| simulated low-rank | 688 | 1.739 | 233 | 208 | 6 | 107 |

Sampler health:

| model | draws | divergences | max treedepth | treedepth hits | min E-BFMI | mean accept |
|---|---:|---:|---:|---:|---:|---:|
| real diagonal | 4,000 | 0 | 8 | 0 | 0.561 | 0.953 |
| simulated diagonal | 4,000 | 0 | 8 | 0 | 0.551 | 0.951 |
| simulated low-rank | 4,000 | 0 | 8 | 0 | 0.586 | 0.951 |

PSIS-LOO was available, but the Pareto-k diagnostics were poor, so the LOO values should not be used as decisive model-selection evidence:

| model | elpd_loo | looic | Pareto k > 0.7 | Pareto k > 1 | max Pareto k |
|---|---:|---:|---:|---:|---:|
| real diagonal | -224,673.4 | 449,346.9 | 4,354 | 3,088 | 3.82 |
| simulated diagonal | -201,763.1 | 403,526.1 | 4,529 | 2,878 | 2.88 |
| simulated low-rank | -196,186.8 | 392,373.7 | 4,215 | 2,104 | 2.08 |

The simulated low-rank model had better nominal predictive performance than the simulated diagonal model on simulated low-rank data, as expected. However, because thousands of observations had Pareto k > 0.7 and many had Pareto k > 1, the PSIS approximation is unstable. A refit-based K-fold comparison would be more defensible.

## Final Interpretation

The real diagonal additive model is acceptable as a descriptive thesis model for variance decomposition. It has clean sampler behavior and only mild convergence warnings, and it gives interpretable player, season, and residual variance shares.

The simulated diagonal fit is a successful sanity check. Even when the data were generated from a low-rank process, the diagonal model recovered player and residual variance diagonals very well and season variance reasonably, though less strongly.

The current low-rank recovery model should not be used for direct loading interpretation. Its sampler did not diverge, but many loading parameters had poor Rhat and very low ESS. Its covariance diagonal results are partly informative, especially for player and residual components, but the loading structure and season-block interpretation are not reliable enough.

For final reporting, use the diagonal additive real-data model as the main descriptive result, report the simulation checks as validation, and present the low-rank model as an attempted richer model that requires stronger identifiability constraints before interpretation.
