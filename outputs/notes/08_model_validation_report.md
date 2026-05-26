# Model Validation Report

Generated: 2026-05-14 13:22:38 CEST

## Artifact Status

- Fit files present: 3 / 3
- Posterior summaries present: 3 / 3

## Posterior Health

# A tibble: 3 × 10
  model_id        status n_parameters max_rhat n_rhat_over_1_01 n_rhat_over_1_05
  <chr>           <chr>         <int>    <dbl>            <int>            <int>
1 03_real_diagon… avail…          531     1.03               13                0
2 05_sim_diagona… avail…          531     1.01                3                0
3 06_sim_lowrank… avail…          688     1.74              233              208
# ℹ 4 more variables: min_ess_bulk <dbl>, min_ess_tail <dbl>,
#   n_ess_bulk_lt_400 <int>, n_ess_tail_lt_400 <int>

## Sampler Health

# A tibble: 3 × 8
  model_id         status n_draws n_divergent max_treedepth n_max_treedepth_hits
  <chr>            <chr>    <int>       <int>         <dbl>                <int>
1 03_real_diagona… avail…    4000           0             8                    0
2 05_sim_diagonal… avail…    4000           0             8                    0
3 06_sim_lowrank_… avail…    4000           0             8                    0
# ℹ 2 more variables: min_ebfmi <dbl>, mean_accept_stat <dbl>

## Simulation Recovery Outputs

- Diagonal misspecification recovery tables are written when `05_sim_diagonal_additive` posterior summaries exist.
- Low-rank loading and variance recovery tables are written when `06_sim_lowrank_recovery` posterior summaries exist.
- LOO summaries are written when fit objects and the `loo` package are available.

Primary CSV outputs:

- `outputs/diagnostics/08_posterior_health_summary.csv`
- `outputs/diagnostics/08_sampler_health_summary.csv`
- `outputs/tables/08_sim_diagonal_recovery_metrics.csv`
- `outputs/tables/08_sim_lowrank_loading_recovery_metrics.csv`
- `outputs/tables/08_sim_lowrank_variance_recovery_metrics.csv`
- `outputs/tables/08_loo_summary.csv`
