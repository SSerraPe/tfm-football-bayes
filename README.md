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

---

<!-- CLAUDE: PIPELINE JOURNAL AUTO-UPDATE INSTRUCTIONS

FILE:    pipeline_journal.qmd  (model root, same directory as run_pipeline.R)
RENDER:  quarto render pipeline_journal.qmd   (run from model root)
OUTPUT:  pipeline_journal.pdf  (appears in model root alongside the .qmd)

TRIGGER: Append a new entry whenever any of the following appear:
  - new files in outputs/diagnostics/  (e.g. *_diagnostic_summary.csv, *_sampler_diagnostics.csv)
  - new files in outputs/tables/       (e.g. *_posterior_summary.csv, *_recovery_metrics.csv)
  - new PNG files in outputs/figures/  (e.g. new chain plots, loading heatmaps)

HOW TO ADD AN ENTRY — step by step:
  1. Read this README note for context.
  2. Read the new output CSV/PNG files to extract the results for the entry.
  3. Open pipeline_journal.qmd and scroll to the VERY END of the file.
  4. Append a new section following the "Entry Template" at the bottom of that file.
     Use the next sequential entry number and today's date.
  5. NEVER delete or edit any previous entry — all entries are permanent journal records.
  6. Run:  quarto render pipeline_journal.qmd
     to produce the updated pipeline_journal.pdf.

ENTRY CONTENT CHECKLIST (fill each subsection from the files listed):
  [ ] Configuration table     — seed, chains, iter_warmup, iter_sampling, adapt_delta, max_treedepth
  [ ] Stages run              — stage numbers + one-line purpose each
  [ ] Sampler health table    — from outputs/diagnostics/*_diagnostic_summary.csv or
                                outputs/diagnostics/11_posterior_health_new_model.csv
  [ ] Recovery metrics table  — from outputs/tables/11_sim_recovery_metrics.csv (sim stages only)
  [ ] Variance decomp table   — from outputs/tables/11_real_variance_decomposition.csv (real stages)
  [ ] Variance decomp plot    — outputs/figures/11_real_variance_decomposition_lowrank.png
  [ ] Loading heatmap         — outputs/figures/11_real_lambda_a_heatmap.png
  [ ] Top loadings plot       — outputs/figures/11_real_lambda_a_top_loadings.png
  [ ] Group radar             — outputs/figures/11_real_lambda_a_group_radar.png
  [ ] Chain trace plots       — outputs/figures/11_chain_plots/{model_tag}/lambda_a_diag_trace.png
  [ ] Chain density plots     — outputs/figures/11_chain_plots/{model_tag}/lambda_a_diag_density.png
  [ ] Bugs fixed              — table of any scripts fixed during this run
  [ ] Notes                   — convergence verdict, football interpretation, next steps

IMPORTANT RULES:
  - Entries are append-only. Never modify or delete Entry 1 or any prior entry.
  - Use the entry template at the bottom of pipeline_journal.qmd as the structural guide.
  - If a file does not exist yet (e.g. stage 12 was not run), simply omit that subsection.
  - Re-render the PDF after every append so pipeline_journal.pdf stays current.
-->

<sub>*`pipeline_journal.qmd` / `pipeline_journal.pdf` — living results journal, auto-updated by Claude on each pipeline run.*</sub>
