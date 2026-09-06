# Hierarchical Bayesian Player Performance Model — La Liga

A Bayesian hierarchical factor model that decomposes La Liga player-season performance into
**player effects** (stable identity), **season effects** (league-wide shifts), and **residual
noise**, fitted with Stan/NUTS via `cmdstanr`. Companion Master's thesis (TFM, MESIO-UPC):
source in `tfm_latex/`.

---

## Mathematical Model

The additive decomposition:

$$\mathbf{y}_n = \mathbf{a}_{i[n]} + \mathbf{b}_{j[n]} + \boldsymbol{\varepsilon}_n$$

where $\mathbf{y}_n \in \mathbb{R}^P$ is a vector of $P$ Z-scaled performance features for
player-season row $n$, $\mathbf{a}_i$ is a player effect, $\mathbf{b}_j$ is a season effect, and
$\boldsymbol{\varepsilon}_n$ is the residual. There is no explicit global mean: features are
Z-scaled across the full sample before fitting, so $\boldsymbol{\mu} \approx \mathbf{0}$.

Covariance structure (asymmetric by design — see `CLAUDE.md` for the full rationale):

| Component | Structure | Notes |
|-----------|-----------|-------|
| $\Sigma_a$ (player) | $\Lambda_a\Lambda_a' + \Psi_a$ | Low-rank ($K$) + diagonal uniqueness; $\Lambda_a$ constrained lower-triangular, positive-diagonal for identification |
| $\Sigma_b$ (season) | $\mathrm{diag}(\sigma_b^2)$ | Diagonal — a low-rank season covariance is under-determined with only $S=12$ seasons (confirmed by preliminary fits: R-hat > 1.70) |
| $\Sigma_e$ (residual) | $\mathrm{diag}(\sigma_e^2)$, Student-$t$ | Heavy-tailed ($\nu$ estimated); optionally scaled by $\sqrt{\bar m/m_n}$ to account for playing-time heteroscedasticity |

Rank $K=3$ is the production choice, selected on qualitative grounds (residual adequacy, common-
variance share and reliable-loader count of the PCs of $\Lambda_a\Lambda_a'$ — see
`tfm_latex/methodology.tex`, §9) rather than by cross-validation, which was found inconclusive
across $K \in \{2,\dots,6\}$.

---

## Data

- **Source:** Wyscout La Liga match-level event data, via a professional Redshift warehouse
- **Sample:** Outfield players only (goalkeepers excluded by match-level activity markers),
  $\geq 600$ minutes played in a season
- **Observations:** $N=4{,}586$ player-seasons across $I=1{,}529$ players and $S=12$ seasons
  (2011/12–2022/23)
- **Features:** $P=48$ Z-scaled per-90 and success-rate statistics (5 log$(1+x)$-transformed,
  18 square-root-transformed, 25 untransformed, before Z-scaling)

---

## Pipeline Stages

The pipeline has grown to more than 30 numbered stages; only the headline ones are listed here —
see `CLAUDE.md` §4 for the complete map (including post-processing, diagnostics, and ablation
stages) and `docs/journal/pipeline_journal.qmd` for the full narrative history.

| Stage | Script | Description |
|-------|--------|-------------|
| 00–02 | `00_build_longitudinal_dataset.R` … `02_prepare_model_objects.R` | Raw data → cleaned, transformed, Z-scaled `Y` + Stan data objects |
| 02b | `02b_descriptive_pca_fa_cmds.R` | Descriptive PCA / factor analysis / classical MDS (feeds thesis Chapter 3) |
| 03 | `03_fit_real_diagonal_additive.R` | Diagonal ($K=0$) baseline |
| 10 | `10_fit_real_lowrank_a_diag_b.R` | Low-rank $\Sigma_a$ + diagonal $\Sigma_b$, Normal residuals ($K=2$) |
| 18 | `18_fit_real_t_errors.R` | Student-$t$ residuals ($K=2$) |
| 28 | `28_fit_k_t_prod.R` | Student-$t$, $K=3$ (production rank) |
| 34 | `34_fit_minutes_scaled_t.R` | Student-$t$, $K=3$, minutes-scaled residual variance (current production model) |
| 11–16 | loading/ICC/interpretation stages | PCA post-processing (rotation-invariant loadings), variance decomposition (ICC), player/season profiles |

---

## Quick Start

```sh
# Data preparation (stages 00-02) — no Stan required
Rscript run_pipeline.R 00 01 02

# Descriptive analysis (no Stan required)
Rscript scripts/02b_descriptive_pca_fa_cmds.R

# Smoke test with fast sampling (2 chains, 200 iterations)
BFA_RUN_STAN=true BFA_CHAINS=2 BFA_ITER_WARMUP=200 BFA_ITER_SAMPLING=200 \
  Rscript run_pipeline.R 10

# Full production fit (K=3, Student-t, minutes-scaled) — several hours
BFA_RUN_STAN=true Rscript run_pipeline.R 34

# Analysis only (reuses existing fits, minutes)
Rscript run_pipeline.R 11 12 13 14 15 16
```

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|--------------|
| `BFA_RUN_STAN` | `false` | Set `true` to run Stan sampling |
| `BFA_REUSE_FIT` | `true` | Reuse existing CSV fits if found |
| `BFA_CHAINS` | `4` | Number of MCMC chains |
| `BFA_PARALLEL_CHAINS` | `4` | Chains running in parallel |
| `BFA_ITER_WARMUP` | `1000` | Warmup iterations per chain |
| `BFA_ITER_SAMPLING` | `1000` | Sampling iterations per chain |
| `BFA_ADAPT_DELTA` | `0.95` | Target acceptance probability |
| `BFA_MAX_TREEDEPTH` | `12` | Max NUTS tree depth |
| `BFA_SIM_RANK_A` | `2` | Player factor rank $K$ used by simulation stages (production fits set this explicitly, e.g. `3` for stage 28/34) |
| `BFA_SIGMA_FLOOR` | `0.05` | Minimum residual/uniqueness scale |
| `BFA_SEED` | `20260513` | Random seed |

---

## Project Structure

```
model/
├── data/
│   ├── raw/          # Original Wyscout CSV (not committed)
│   ├── interim/      # Intermediate aggregations (not committed)
│   └── processed/    # Final scaled data and model objects
├── src/               # Shared helper functions (bootstrap.R, config.R, stan_helpers.R, ...)
├── stan/               # Stan model files, incl. stan/prototypes/ (early, superseded drafts)
├── stan_data/          # Serialised Stan data objects
├── fits/               # Lightweight RDS fit pointers (CSV draws under fits/csv/, gitignored)
├── scripts/            # Numbered pipeline stages (00 through 37+)
├── outputs/
│   ├── tables/         # CSV result tables
│   ├── figures/        # PNG plots
│   ├── diagnostics/    # Sampler and posterior health
│   └── notes/          # Generated stage notes (script output)
├── docs/                # All written documentation — journal, professor summary, thesis,
│                        #   model docs, ideas, analysis notes, runtime working doc
├── archive/              # Superseded/historical material, incl. the migrated exploratory
│                        #   p=23 pipeline and Bayesian-MDS prototype (see its own README)
├── rank comparison/     # Standalone CV-ELPD rank-selection methodology
├── notebooks/           # R Markdown exploratory notebooks
├── tfm_latex/            # The thesis document itself (LaTeX, UPC/FME template)
└── run_pipeline.R       # Pipeline entry point
```

---

## Key Results

| Component | Finding |
|-----------|---------|
| Low-rank vs. diagonal $\Sigma_a$ | Strongly favoured (ΔELPD ≈ +3,724, PSIS-LOO, pre-rebuild data) |
| Student-$t$ vs. Normal residuals | Strongly favoured (ΔELPD ≈ +11,668); $\hat\nu$ moved 2.95 → 4.9 after GK exclusion/transforms → 5.7–6.0 after minutes-scaling |
| Rank $K$ | $K^*=3$, on qualitative grounds (see Mathematical Model above) |
| Minutes heteroscedasticity | Confirmed and corrected — low-minute share of worst residuals dropped from 48.9% to 8.3% after minutes-scaling |
| Player archetypes | 3 canonical factor directions (PCA of $\Lambda_a\Lambda_a'$), interpretable as defensive-work-rate vs. attacking output, technical/pass vs. aerial/physical, and possession quality vs. wide/direct play |

See `docs/professor/professor_summary.pdf` for the full results write-up and `CLAUDE.md` §8 for
the detailed, currently-open items.

---

## Reproducibility Note

Stan posterior CSV draws (tens of GB) are excluded from version control due to size.
Lightweight fit pointers in `fits/*.rds` allow reconstruction via `load_cmdstan_fit()`
(`src/stan_helpers.R`), which falls back to path discovery if the stored paths are stale.

---

## Requirements

- R ≥ 4.2
- CmdStan ≥ 2.32 (install via `cmdstanr::install_cmdstan()`)
- R packages: `cmdstanr`, `posterior`, `bayesplot`, `ggplot2`, `dplyr`, `tidyr`, `readr`, `tibble`,
  `scales`, `psych` (factor analysis), `ggrepel` (label annotations)
- Optional: `quarto` (PDF reports), `bayMDS` (only for the archived Bayesian-MDS prototype)

---

<!-- CLAUDE: PIPELINE JOURNAL AUTO-UPDATE INSTRUCTIONS

FILE:    docs/journal/pipeline_journal.qmd  (moved here 2026-06-26)
RENDER:  quarto render docs/journal/pipeline_journal.qmd
OUTPUT:  docs/journal/pipeline_journal.pdf  (gitignored; regenerated on demand)
NOTE:    the .qmd sets its knit root two levels up to the project root, so all
         outputs/... paths resolve when rendered from anywhere.

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
  6. Run:  quarto render docs/journal/pipeline_journal.qmd
     to produce the updated docs/journal/pipeline_journal.pdf.

ENTRY CONTENT CHECKLIST (fill each subsection from the files listed, when applicable):
  [ ] Configuration table     — seed, chains, iter_warmup, iter_sampling, adapt_delta, max_treedepth
  [ ] Stages run              — stage numbers + one-line purpose each
  [ ] Sampler health table    — from outputs/diagnostics/*_diagnostic_summary.csv
  [ ] Recovery metrics table  — from outputs/tables/*_recovery_metrics.csv (sim stages only)
  [ ] Variance decomp table   — from outputs/tables/*_variance_decomposition.csv (real stages)
  [ ] Key figures              — reference the relevant outputs/figures/*.png
  [ ] Bugs fixed              — table of any scripts fixed during this run
  [ ] Notes                   — convergence verdict, football interpretation, next steps

IMPORTANT RULES:
  - Entries are append-only. Never modify or delete any prior entry.
  - Use the entry template at the bottom of pipeline_journal.qmd as the structural guide.
  - If a file does not exist yet, simply omit that subsection.
  - Re-render the PDF after every append so pipeline_journal.pdf stays current.
-->

<sub>*`pipeline_journal.qmd` / `pipeline_journal.pdf` — living results journal, auto-updated by Claude on each pipeline run.*</sub>
