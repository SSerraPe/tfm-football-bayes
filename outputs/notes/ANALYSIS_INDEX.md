# Analysis Output Index

Quick reference for all outputs in `outputs/`.

---

## Tables (`outputs/tables/`)

### Data pipeline (stages 00-02)
| File | Description |
|------|-------------|
| `00_raw_schema_summary.csv` | Column types and missingness in raw Wyscout CSV |
| `00_aggregation_rules.csv` | Aggregation logic applied per variable |
| `00_aggregation_summary.csv` | Row counts before/after match-to-season aggregation |
| `00_player_counts_by_season.csv` | Player count per season after aggregation |
| `01_cleaning_summary.csv` | Rows dropped at each filtering step |
| `01_feature_missingness_summary.csv` | Missingness rate per feature before/after imputation |
| `01_feature_variance_summary.csv` | Variance check per feature |
| `01_feature_inclusion_exclusion_audit.csv` | Final feature inclusion decision |
| `01_selected_variables.csv` | The 53 selected feature names |
| `02_feature_matrix_audit_summary.csv` | Rank, max correlation of scaled feature matrix |
| `02_feature_high_correlation_pairs.csv` | Feature pairs with correlation > 0.95 |
| `02_model_object_summary.csv` | N, P, I, S dimensions |
| `02_player_counts_by_season.csv` | Player counts per season in model |
| `02_player_recurrence_diagnostics.csv` | How many seasons each player appears |

### Posterior summaries (stages 03, 05, 06, 09, 10)
| File | Description |
|------|-------------|
| `03_real_diagonal_additive_posterior_summary.csv` | Real data diagonal model: 531 parameters |
| `05_sim_diagonal_additive_posterior_summary.csv` | Simulated data diagonal model |
| `06_sim_lowrank_recovery_posterior_summary.csv` | Full low-rank recovery model (FAILED: Rhat 1.739) |
| `09_sim_lowrank_a_diag_b_posterior_summary.csv` | New model on simulated data |
| `10_real_lowrank_a_diag_b_posterior_summary.csv` | New model on real La Liga data |

### Simulation truth (stage 04)
| File | Description |
|------|-------------|
| `04_simulation_summary.csv` | Simulation parameters and dimensions |
| `04_simulation_truth_loadings.csv` | Known Lambda_a and Lambda_b used to generate data |
| `04_simulation_truth_variance_components.csv` | Known Sigma_a_diag, Sigma_b_diag, Sigma_e_diag |

### Variance decomposition (stages 07, 11)
| File | Description |
|------|-------------|
| `07_real_variance_decomposition.csv` | Player/season/residual proportions — real diagonal model |
| `07_real_posterior_variance_components.csv` | Posterior mean variances per feature |
| `07_real_empirical_variance_decomposition.csv` | Empirical decomposition (data-based, no model) |
| `11_real_variance_decomposition.csv` | Player/season/residual proportions — new model |

### Recovery metrics (stages 07, 08, 11)
| File | Description |
|------|-------------|
| `07_sim_diagonal_variance_vs_truth.csv` | Diagonal model variance recovery vs truth |
| `07_sim_lowrank_recovery_loadings.csv` | Loading estimates vs truth (low-rank model) |
| `07_sim_lowrank_recovery_variances.csv` | Variance recovery for low-rank model |
| `08_sim_diagonal_recovery_metrics.csv` | Aggregated MAE/RMSE/correlation for diagonal recovery |
| `08_sim_lowrank_loading_recovery_metrics.csv` | Loading recovery metrics (failed model) |
| `08_sim_lowrank_variance_recovery_metrics.csv` | Variance recovery metrics |
| `08_sim_diagonal_recovery_by_feature.csv` | Feature-level diagonal recovery |
| `08_sim_lowrank_loading_recovery_by_parameter.csv` | Parameter-level loading recovery |
| `08_sim_lowrank_variance_recovery_by_feature.csv` | Feature-level variance recovery |
| `11_sim_recovery_metrics.csv` | Recovery metrics for new model vs simulation truth |

### Loading tables (stage 11)
| File | Description |
|------|-------------|
| `11_sim_lambda_a_loading_table.csv` | Lambda_a loadings with CI — simulated data fit |
| `11_real_lambda_a_loading_table.csv` | Lambda_a loadings with CI — real La Liga fit |

### Cross-validation and model comparison (stages 08, 20)
| File | Description |
|------|-------------|
| `08_loo_summary.csv` | PSIS-LOO results (UNRELIABLE: Pareto k too high) |
| `20_loo_diagonal_vs_lowrank.csv` | LOO comparison: diagonal (stage 03) vs low-rank (stage 10) |
| `20_icc_diagonal_vs_lowrank.csv` | ICC_player per feature, both models side by side |

### Posterior predictive checks (stage 19)
| File | Description |
|------|-------------|
| `19_correlation_ppc.csv` | Observed vs model-implied correlation per feature pair |

### Rank selection (stage 17)
| File | Description |
|------|-------------|
| `outputs/rank_selection/real_data_elpd.csv` | Held-out ELPD by rank (K=0..4) |
| `outputs/rank_selection/real_data_holdouts.csv` | Fold structure and held-out counts |

### Football interpretation (stage 12)
| File | Description |
|------|-------------|
| `12_player_factor_scores.csv` | Posterior mean factor scores per player |

---

## Figures (`outputs/figures/`)

### Variance decomposition
| File | Description |
|------|-------------|
| `07_real_empirical_variance_decomposition.png` | Stacked bar: player/season/residual shares (real data) |
| `11_real_variance_decomposition_lowrank.png` | Same for new low-rank Sigma_a model |

### Recovery (simulation)
| File | Description |
|------|-------------|
| `07_sim_diagonal_covariance_diag_vs_truth.png` | Diagonal variance recovery heatmap |
| `07_sim_lowrank_loadings_vs_truth.png` | Loading estimates vs truth scatter |
| `07_sim_lowrank_variances_vs_truth.png` | Variance component recovery |
| `08_sim_diagonal_recovery.png` | Aggregated recovery metrics |
| `08_sim_lowrank_loading_recovery.png` | Loading-level recovery heatmap |
| `11_sim_lambda_a_vs_truth.png` | New model: Lambda_a posterior mean vs truth |
| `11_sim_recovery_metrics.png` | Recovery metrics summary for new model |

### Loading visualizations (stage 11)
| File | Description |
|------|-------------|
| `11_sim_lambda_a_heatmap.png` | Features × factors heatmap (simulated data) |
| `11_real_lambda_a_heatmap.png` | Features × factors heatmap (real La Liga data) |
| `11_sim_lambda_a_top_loadings.png` | Top loadings bar chart with sign (simulated) |
| `11_real_lambda_a_top_loadings.png` | Top loadings bar chart with sign (real data) |
| `11_real_lambda_a_group_radar.png` | Mean |loading| by feature group per factor |

### Football interpretation (stage 12)
| File | Description |
|------|-------------|
| `12_real_factor_group_radar.png` | Group-level radar for each factor |
| `12_real_loading_comparison.png` | Factor loadings per feature, side by side |
| `12_player_archetype_scatter.png` | Player positions in factor score space |

### Posterior health (stage 08)
| File | Description |
|------|-------------|
| `08_posterior_health_summary.png` | Rhat and ESS summary across all models |

### Posterior predictive checks (stage 19)
| File | Description |
|------|-------------|
| `19_correlation_ppc_scatter.png` | Observed vs model-implied correlation — scatter plot |
| `19_residual_distribution_check.png` | Residual histograms vs Normal(0, σ_e) for selected features |

### Model comparison (stage 20)
| File | Description |
|------|-------------|
| `20_icc_diagonal_vs_lowrank.png` | ICC_player scatter: diagonal (stage 03) vs low-rank (stage 10) |

### Chain plots (stage 11)
| Directory | Description |
|-----------|-------------|
| `11_chain_plots/09_sim/` | Trace and density plots — simulated data fit |
| `11_chain_plots/10_real/` | Trace and density plots — real data fit |

---

## Diagnostics (`outputs/diagnostics/`)

| File | Description |
|------|-------------|
| `03/05/06_*_diagnostic_summary.csv` | Quick per-model health summary |
| `03/05/06_*_sampler_diagnostics.csv` | Per-draw sampler diagnostics |
| `08_posterior_health_summary.csv` | Rhat, ESS across all models |
| `08_sampler_health_summary.csv` | Divergences, E-BFMI, acceptance rates |
| `11_posterior_health_new_model.csv` | Health summary for stages 09 and 10 |

---

## Notes (`outputs/notes/`)

| File | Description |
|------|-------------|
| `08_model_validation_report.md` | Artifact availability and health status |
| `09_validation_interpretation.md` | What results mean for the thesis |
| `10_comprehensive_validation_analysis.md` | Full technical validation (1640 lines) |
| `12_factor_interpretation.md` | Football interpretation of factor loadings |
| `note_mu_centering.md` | Why μ is redundant (data already Z-scaled) |
| `note_explained_variability.md` | ICC_player and ICC_season interpretation |
| `note_player_interpretation.md` | PCA of ΛΛ', factor meanings, player clustering |
| `note_rank_selection_methodology.md` | CV-ELPD rank selection: how and why |
| `note_rank_selection_results.md` | Selected rank K* (written after stage 17 runs) |
| `note_t_errors.md` | Student-t residuals: motivation and interpretation |
| `note_posterior_predictive_checks.md` | What the PPCs test and how to read them |
| `note_diagonal_vs_lowrank.md` | LOO comparison result: low-rank wins by Δ=3724 (19.9 SE) |
| `ANALYSIS_INDEX.md` | This file |
