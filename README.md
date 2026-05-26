# Hierarchical Bayesian Player Performance Model — La Liga

A complete Bayesian modelling pipeline that decomposes La Liga player-season performance into **player effects**, **season effects**, and **residual noise** using additive hierarchical models fitted with Stan/NUTS.

---

## Mathematical Model

The core additive decomposition:

$$\mathbf{y}_n = \boldsymbol{\mu} + \mathbf{a}_{i[n]} + \mathbf{b}_{j[n]} + \boldsymbol{\varepsilon}_n$$

where $\mathbf{y}_n \in \mathbb{R}^P$ is a vector of P scaled performance features for player-season row $n$, $\mathbf{a}_i$ is a player effect, $\mathbf{b}_j$ is a season effect, and $\boldsymbol{\varepsilon}_n$ is residual noise.

Four covariance specifications are implemented:

| Model | $\Sigma_a$ (player) | $\Sigma_b$ (season) | Stage |
|-------|-------------------|-------------------|-------|
| Diagonal baseline | $\text{Diag}(\sigma_a^2)$ | $\text{Diag}(\sigma_b^2)$ | 03 |
| Sim diagonal (misspecification check) | $\text{Diag}(\sigma_a^2)$ | $\text{Diag}(\sigma_b^2)$ | 05 |
| Full low-rank recovery (failed) | $\Lambda_a\Lambda_a' + \Psi_a$ | $\Lambda_b\Lambda_b' + \Psi_b$ | 06 |
| **Low-rank $\Sigma_a$ + diagonal $\Sigma_b$** | $\Lambda_a\Lambda_a' + \Psi_a$ | $\text{Diag}(\sigma_b^2)$ | **09/10** |

The new model (stages 09–10) is the key contribution: it captures player-level latent factors (archetypes) while keeping season effects simple — the right tradeoff given only 12 seasons.

---

## Data

- **Source:** Wyscout La Liga player-match event data
- **Sample:** Non-goalkeeper players with ≥450 minutes in a season
- **Observations:** 4,944 player-season rows across 1,650 players and 12 seasons
- **Features:** 53 scaled (Z-scored) per-90 and success-rate statistics

---

## Pipeline Stages

| Stage | Script | Description |
|-------|--------|-------------|
| 00 | `00_build_longitudinal_dataset.R` | Aggregate match rows to player-season panel |
| 01 | `01_read_clean_feature_engineering.R` | Clean, filter, engineer 53 features |
| 02 | `02_prepare_model_objects.R` | Scale, build indices and Stan data |
| 03 | `03_fit_real_diagonal_additive.R` | Fit diagonal model to real data (main result) |
| 04 | `04_simulate_lowrank_additive_data.R` | Simulate data with known low-rank truth |
| 05 | `05_fit_sim_diagonal_additive.R` | Diagonal model on simulated data (misspecification check) |
| 06 | `06_fit_sim_lowrank_recovery.R` | Full low-rank recovery (failed — poor season mixing) |
| 07 | `07_diagnostics_and_recovery_plots.R` | Recovery tables and plots |
| 08 | `08_validate_all_models.R` | Comprehensive validation (sampler, posterior, LOO) |
| **09** | `09_fit_sim_lowrank_a_diag_b.R` | **New model on simulated data (validation)** |
| **10** | `10_fit_real_lowrank_a_diag_b.R` | **New model on real data (factor discovery)** |
| **11** | `11_loading_analysis_and_chains.R` | **Loading tables, chain plots, recovery metrics** |
| **12** | `12_football_interpretation.R` | **Factor interpretation, player archetypes** |

---

## Quick Start

```sh
# Data preparation (stages 00-02) — no Stan required
Rscript model/run_pipeline.R 00 01 02

# Smoke test with fast sampling (2 chains, 200 iterations)
BFA_RUN_STAN=true BFA_CHAINS=2 BFA_ITER_WARMUP=200 BFA_ITER_SAMPLING=200 \
  Rscript model/run_pipeline.R 09

# Full new-model pipeline (requires prior fits from 03-08)
BFA_RUN_STAN=true Rscript model/run_pipeline.R 09 10 11 12

# Analysis only (if fits already exist)
Rscript model/run_pipeline.R 11 12
```

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BFA_RUN_STAN` | `false` | Set `true` to run Stan sampling |
| `BFA_REUSE_FIT` | `true` | Reuse existing CSV fits if found |
| `BFA_CHAINS` | `4` | Number of MCMC chains |
| `BFA_PARALLEL_CHAINS` | `4` | Chains running in parallel |
| `BFA_ITER_WARMUP` | `1000` | Warmup iterations per chain |
| `BFA_ITER_SAMPLING` | `1000` | Sampling iterations per chain |
| `BFA_ADAPT_DELTA` | `0.95` | Target acceptance probability |
| `BFA_MAX_TREEDEPTH` | `12` | Max NUTS tree depth |
| `BFA_SIM_RANK_A` | `2` | Player factor rank Q_a |
| `BFA_SIM_RANK_B` | `1` | Season factor rank (stages 04/06 only) |
| `BFA_SEED` | `20260513` | Random seed |

---

## Project Structure

```
model/
├── data/
│   ├── raw/          # Original Wyscout CSV (not committed)
│   ├── interim/      # Intermediate aggregations (not committed)
│   └── processed/    # Final scaled data and model objects
├── src/              # Shared helper functions
│   ├── bootstrap.R
│   ├── config.R
│   ├── stan_helpers.R
│   ├── diagnostics_helpers.R
│   ├── loading_visualization.R
│   └── ...
├── stan/             # Stan model files
│   ├── additive_diagonal.stan
│   ├── additive_lowrank_recovery.stan
│   └── additive_lowrank_a_diag_b.stan
├── stan_data/        # Serialised Stan data objects
├── fits/             # Lightweight RDS fit pointers
├── scripts/          # Numbered pipeline stages 00-12
├── outputs/
│   ├── tables/       # CSV result tables
│   ├── figures/      # PNG plots
│   ├── diagnostics/  # Sampler and posterior health
│   └── notes/        # Markdown documents
├── reports/          # Compiled thesis document (Quarto)
├── docs/             # Technical model documentation
├── notebooks/        # R Markdown exploratory notebooks
└── run_pipeline.R    # Pipeline entry point
```

---

## Key Results

| Component | Finding |
|-----------|---------|
| Player effects | Dominant (60–88%) for per-90 volume features |
| Season effects | Moderate (15–38%) for tactical/style features |
| Residual | Dominant (80–98%) for success rates |
| New model (09/10) | Low-rank player factors with diagonal season effects |

---

## Reproducibility Note

Stan posterior CSV draws (~12 GB) are excluded from version control due to size.  
Lightweight fit pointers in `fits/*.rds` allow reconstruction via `cmdstanr::as_cmdstan_fit()`.

---

## Requirements

- R ≥ 4.2
- CmdStan ≥ 2.32 (install via `cmdstanr::install_cmdstan()`)
- R packages: `cmdstanr`, `posterior`, `bayesplot`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `tibble`, `scales`
- Optional: `ggrepel` (player label annotations), `quarto` (PDF report)
