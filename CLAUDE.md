# CLAUDE.md — TFM Football Analytics Model

Durable context for future sessions. Read this first, then check the latest entry in
`pipeline_journal.qmd` for where work was last left off. Companion document:
`PROJECT_RUNDOWN_AND_RUNTIME.txt` (plain-language walkthrough + runtime working doc).

---

## 1. Project goal

Master's thesis (TFM, MESIO UPC). Build a **hierarchical Bayesian factor model** of La Liga
player-season performance that (a) decomposes each performance metric into stable player
identity, season-wide effects, and noise, and (b) recovers interpretable latent "player
archetype" dimensions from the data. Supervised by a Bayesian-statistics professor; modelling
choices are driven by his feedback (recorded in the journal, Entries 1–4).

The user is a professional football analyst, so the latent factors are expected to map onto
recognisable football archetypes, and classical/frequentist factor analysis in that domain
typically uses ~10 factors (relevant to the rank-selection discussion).

---

## 2. Model architecture (current)

Additive decomposition of the P-dimensional observation for player-season `n`:

```
y_n = a_{i[n]} + b_{j[n]} + ε_n          n = 1..N
```

- `a_i ∈ R^P` — player effect (stable identity), `i = 1..I`
- `b_j ∈ R^P` — season effect, `j = 1..S`
- `ε_n ∈ R^P` — residual

Dimensions: **N = 4944** player-seasons, **I = 1650** players, **S = 12** seasons (2011–2023),
**P = 53** Z-scaled performance features.

### Covariance structure

```
Σ_a = Λ_a Λ_a' + Ψ_a       (low-rank factor part + diagonal "uniqueness")
Σ_b = diag(σ_b²)           (diagonal season covariance)
Σ_e = diag(σ_e²)           (diagonal residual covariance)
```

- `Λ_a ∈ R^{P×K}` — factor loadings; production fits use **K = 2** (a proof-of-concept rank,
  not the final choice — see rank selection below).
- `Ψ_a = diag(ψ_a²)` — per-feature uniqueness not explained by the K shared factors.
- Player effect built non-centred: `a_i = Λ_a η_i + ψ_a ⊙ z_{u,i}`, with `η_i` (factor scores)
  and `z_{u,i}` (uniqueness scores) standard normal.

### Likelihood

- **Normal residuals** (stage 10): `ε_np ~ N(0, σ_e,p)`.
- **Student-t residuals** (stage 18, preferred): `ε_np ~ t(ν, 0, σ_e,p)`, `ν ~ Gamma(2, 0.1)`.
  Fitted `ν ≈ 2.95` (95% CI 2.91–2.99) → very heavy tails.

### Identification

- `Λ_a` is **lower-triangular with positive diagonal** (`make_lower_tri_loadings` in the Stan
  files; `make_plt()` in `src/stan_helpers.R`). This resolves the rotational non-identifiability
  of the factor model. Loadings are therefore **not directly interpretable** — interpretation
  comes from a PCA re-orientation of `Λ_a Λ_a'` per draw (stage 13).
- Scores are centred (sum-to-zero across players/seasons) inside the Stan `transformed
  parameters` block; centring is applied to the raw scores, not the effects, to avoid
  `inf − inf = NaN` (see comment block in `stan/additive_lowrank_a_diag_b_t.stan`).

### Variance decomposition / ICC

Per feature `p`, with `Σ_a,pp = Λ_p·Λ_p' + ψ_a,p²`, `Σ_b,pp = σ_b,p²`, `Σ_e,pp = σ_e,p²`:

```
ICC_player   = Σ_a,pp / (Σ_a,pp + Σ_b,pp + Σ_e,pp)
ICC_season   = Σ_b,pp / (...)
ICC_residual = Σ_e,pp / (...)
```

For the t-model, `Var(ε) = σ_e² · ν/(ν−2)` is used (see `generated quantities`).

---

## 3. Key modelling decisions (and why)

| Decision | Rationale |
|---|---|
| Low-rank `Σ_a` (vs diagonal) | LOO-CV strongly favours low-rank: ΔELPD ≈ +3,724 (≈20 SE) over the K=0 diagonal baseline. Shared factor structure is predictively real. |
| Lower-triangular `Λ_a` | Removes rotational non-identifiability; confirmed needed by the professor and by simulation recovery (Entry 3). |
| PCA post-processing of `Λ_aΛ_a'` | The triangular `Λ` is an identification device, not interpretable. Eigendecomposing `ΛΛ'` per draw gives a canonical, rotation-invariant set of loadings + player scores (stage 13). |
| PCA-based initialisation | `build_pca_init()` warm-starts chains from a classical factor-analysis solution of the player means. Materially eases convergence (professor's `make_crossed_init`, adapted). |
| Data Z-scaled, no global `μ` | Features are standardised within season in stage 02, so `μ ≈ 0`. Newer models (stage 18+) drop the `μ` parameter entirely. |
| Student-t residuals | A few player-seasons (breakouts, injury returns, odd tactical roles) are extreme; `ν ≈ 3` and ΔELPD ≈ +11,668 confirm heavy tails. **Default for future fits.** |
| ICC framing of `prop_a/b/e` | Professor reframed the variance proportions as intraclass correlations — the natural "explained variability" quantity for the thesis. |

---

## 4. Pipeline map (`run_pipeline.R` → `scripts/`)

Stages run as isolated `Rscript` processes. `run_pipeline.R 10 11 12` runs a subset in order.

| Stage | Script | Purpose |
|---|---|---|
| 00 | `00_build_longitudinal_dataset.R` | Wyscout raw → longitudinal player-season panel |
| 01 | `01_read_clean_feature_engineering.R` | Cleaning + per90 / rate feature engineering |
| 02 | `02_prepare_model_objects.R` | Build Z-scaled `Y`, indices, `model_objects.rds`, maps |
| 03 | `03_fit_real_diagonal_additive.R` | Diagonal baseline (K=0) on real data |
| 04 | `04_simulate_lowrank_additive_data.R` | Simulate known low-rank truth |
| 05 | `05_fit_sim_diagonal_additive.R` | Diagonal fit on simulated data |
| 06 | `06_fit_sim_lowrank_recovery.R` | Low-rank recovery check on simulated data |
| 07 | `07_diagnostics_and_recovery_plots.R` | Recovery + diagnostics plots |
| 08 | `08_validate_all_models.R` | Cross-model validation |
| 09 | `09_fit_sim_lowrank_a_diag_b.R` | Fit low-rank `Σ_a` + diag `Σ_b` to simulated data |
| 10 | `10_fit_real_lowrank_a_diag_b.R` | **Real-data production fit (Normal, K=2)** |
| 11 | `11_loading_analysis_and_chains.R` | Loadings table, variance decomp, chain diagnostics |
| 12 | `12_football_interpretation.R` | Factor interpretation, radar (Fig 4), archetype scatter |
| 13 | `13_pca_postprocessing.R` | PCA of `ΛΛ'` per draw → loadings + player scores |
| 14 | `14_player_profiles.R` | Player heatmap + PCA scatter |
| 15 | `15_season_profiles.R` | Season effects `b_j` over time (Fig 7) |
| 16 | `16_icc_analysis.R` | ICC decomposition per feature |
| 17 | `17_rank_selection_real.R` | Held-out CV-ELPD for K=0..4 (uses `rank comparison/`) |
| 18 | `18_fit_real_t_errors.R` | **Student-t residual fit (K=2), LOO vs stage 10** |
| 18b | `18_postprocess.R` | LOO + ν post-processing for stage 18 |
| 19 | `19_posterior_predictive_checks.R` | Correlation PPC + residual distribution check |
| 20 | `20_diagonal_comparison.R` | LOO: low-rank vs diagonal `Σ_a` |
| 21 | `21_fit_k4_t.R` | **Untracked/new:** K=4 t-errors fit (not yet in `run_pipeline.R`) |

Note: stage 21 exists as a script but is **not** registered in `run_pipeline.R`'s `stage_map`.

---

## 5. Environment variables (`src/config.R`)

Fitting is gated and configurable via `BFA_*` env vars. Key ones:

- `BFA_RUN_STAN` (default `false`) — must be `true` to actually sample; otherwise stages skip.
- `BFA_REUSE_FIT` (default `true`) — reuse existing CSV fits instead of resampling.
- `BFA_CHAINS` / `BFA_PARALLEL_CHAINS` (4 / 4), `BFA_ITER_WARMUP` / `BFA_ITER_SAMPLING` (1000 / 1000).
- `BFA_ADAPT_DELTA` (0.95), `BFA_MAX_TREEDEPTH` (12), `BFA_SIM_RANK_A` (2 — the K used).
- `BFA_SEED` (20260513), `BFA_SIGMA_FLOOR` (0.05).

Typical re-fit (background, hours): `BFA_RUN_STAN=true BFA_REUSE_FIT=false Rscript run_pipeline.R 10`.
Re-run analysis only (reuses fits, minutes): `Rscript run_pipeline.R 11 12 13 14 15 16`.

---

## 6. Folder / file structure

```
model/
├─ run_pipeline.R            # stage runner (00..20)
├─ scripts/                  # numbered pipeline stages 00..21
├─ src/                      # shared helpers, sourced via bootstrap.R
│  ├─ config.R               # paths, BFA_* env vars, curated feature list
│  ├─ stan_helpers.R         # stan_data builders, build_pca_init(), fit load/persist
│  ├─ loading_visualization.R# feature_group_lookup(), group colours, loading plots
│  ├─ data_helpers.R / feature_helpers.R / diagnostics_helpers.R / simulation_helpers.R
│  ├─ large_fit_helpers.R    # extract_stan_csv_params() — column-pass for 12GB stage-18 CSVs
│  ├─ extract_stan_csv_params.py  # Python column-pass script (called by large_fit_helpers.R)
│  └─ bootstrap.R            # sources config + helpers, check_packages()
├─ stan/                     # Stan model files
│  ├─ additive_diagonal.stan
│  ├─ additive_lowrank_recovery.stan
│  ├─ additive_lowrank_a_diag_b.stan        # main Normal model
│  └─ additive_lowrank_a_diag_b_t.stan      # Student-t model (no mu)
├─ data/processed/           # Y_scaled.csv, model_objects.rds, player_map.csv, season_map.csv
├─ stan_data/                # cached *_stan_data.rds inputs
├─ fits/                     # lightweight fit pointers (*.rds) + fits/csv/<id>/<timestamp>/
├─ outputs/{tables,figures,diagnostics,notes,rank_selection}/
├─ logs/                     # background-run logs (timing lives here)
├─ docs/                     # ← ALL written documentation lives here (consolidated 2026-06-26)
│  ├─ PROJECT_RUNDOWN_AND_RUNTIME.txt
│  ├─ journal/pipeline_journal.qmd        # detailed journal — Entries 1..4
│  ├─ professor/professor_summary.qmd + .pdf   # results summary (shareable PDF kept)
│  ├─ thesis/thesis_technical_document.qmd + .pdf + references.bib
│  ├─ model/   additive_model.{md,tex}, modelling_steps_03_to_08.md
│  ├─ ideas/   ideas.Rmd, TODO.Rmd
│  └─ analysis_notes/  ANALYSIS_INDEX.md, note_*.md, validation writeups
├─ archive/                  # legacy/historical material (not active)
│  ├─ rank_comparison_previous/   # historical predecessor pipeline
│  └─ notes_latex_builds/         # regenerable .tex/.pdf writeups
├─ rank comparison/          # standalone CV-ELPD rank-selection methodology + Stan
└─ CLAUDE.md, README.md      # context anchors (kept at root by convention)
```

Note: only `professor_summary.pdf` and `thesis_technical_document.pdf` are tracked PDFs
(shareable deliverables). All other rendered PDFs are gitignored and regenerated on demand.
The three render-source `.qmd` files set their knit root two levels up to the project root, so
their `outputs/...` includes still resolve after the move.

Fits are stored as **lightweight pointers**: the `.rds` in `fits/` references CmdStan CSV files
under `fits/csv/<model_id>/<timestamp>/`. `load_cmdstan_fit()` (`src/stan_helpers.R`) reconstructs
the fit, with fallback discovery if stored paths are stale (e.g. after moving the project).

---

## 7. Supporting documentation (where to look)

- **`docs/journal/pipeline_journal.qmd`** — the canonical narrative log. Entry 1 (stages 03–12
  baseline → low-rank), Entry 2 (professor feedback: PCA post-proc, profiles, ICC, init), Entry 3
  (LT identification confirmed; stage 10 re-fit), Entry 4 (rank selection, t-errors, PPC, model
  comparison), Entry 5 (diagnostics, screening, transforms, runtime hardening), Entry 6 (GK
  discovery, residual calibration, Lambda_a mixing check, feature decision). **Always read the latest entry when resuming.**
- **`docs/professor/professor_summary.qmd`** — polished results doc; section per professor TODO.
  Figure numbering used by the professor maps to: Fig 1 = ICC scatter, Fig 2 = ELPD-rank,
  Fig 3 = loading heatmap, **Fig 4 = group radar (`12_real_factor_group_radar.png`)**,
  Fig 5 = player scatter, Fig 6 = player heatmap,
  **Fig 7 = season profiles grid (`15_season_profiles_grid.png`)**.
- **`docs/thesis/thesis_technical_document.qmd`** — the formal thesis write-up.
- **`docs/analysis_notes/`** — hand-authored interpretation notes (`ANALYSIS_INDEX.md`,
  `note_*.md`). Generated stage notes stay in `outputs/notes/`.
- **`rank comparison/`** — self-contained rank-selection methodology (`methodology_report.Rmd`,
  `rank_selection.R`, `run_on_real_data.R`, Stan models). Folder name kept (with space) because
  `scripts/17_rank_selection_real.R` references it.
- **`docs/PROJECT_RUNDOWN_AND_RUNTIME.txt`** — plain-language process walkthrough + the runtime
  working doc + per-issue fix notes + Claude-Project upload list.

---

## 8. Current state & known issues (as of 2026-06-30)

**State:** production fits are stage 10 (Normal, K=2) and stage 18 (Student-t, K=2; actual
sampler settings: 500 warmup + 1000 sampling, delta=0.9, max_depth=10). Rank selection
(stage 17) shows CV-ELPD rising monotonically through K=4 with no elbow. A K=4 t-fit
script (stage 21) exists but the fit has not been finalised/registered.

**Implemented Session 5 (2026-06-30):**
- `log_lik` gating: both Stan models gate the N-loop behind `BFA_COMPUTE_LOG_LIK` (default 0).
  Set to 1 for any fit where LOO will be run afterwards.
- Pathfinder mode added to rank-selection: `BFA_RANK_METHOD=pathfinder BFA_RANK_MAX=8` runs
  fast approximate screening for K=5–8 before committing to a full NUTS sweep.
- Stage-18 diagnostics: Python-extracted draws (587 small params, col-pass over 12GB CSV) +
  `posterior::summarise_draws()`; written to `18_posterior_health.csv`. Key finding: Lambda_a
  loadings have ESS ≈ 19 (worst) and Rhat up to 1.185. nu and sigma_e are well-converged (ESS
  1450/1550, Rhat < 1.003). 167 params with Rhat > 1.01 (vs 156 in stage-10).
- Feature analysis: stage 22 (residuals + ICC screen), stage 23 (transform recommendations),
  stage 24 (Frobenius, K-fold ELPD, ICC shift, Pareto-k reliability note).

**Implemented Session 6 (2026-06-30, this session):**
- GK player exclusion confirmed missing from pipeline (only GK-specific *column* names were
  filtered in stage 01; GK player *rows* remained). Confirmed 358 GK rows / 121 GK players
  via `total_gk_saves > 0` join from raw Wyscout CSV. Fix goes in stage 01's eligibility step.
- Large-fit extraction utility: `src/extract_stan_csv_params.py` (Python column-pass) +
  `src/large_fit_helpers.R` (R wrapper `extract_stan_csv_params()`). Resolves the 45-min
  cmdstanr bottleneck for the 12GB stage-18 CSVs going forward.
- Lambda_a Lambda_a' mixing check (Task S): ESCALATED with partial rotation-artifact improvement.
  Diagonal ESS=104, Rhat=1.056; off-diagonal ESS=24, Rhat=1.130. 5.5× improvement vs. raw
  Lambda_a (rotation-labeling contribution) but still below thresholds. GK bimodal structure is
  the leading geometry hypothesis. Task O blocked pending Task K.
- Feature keep/drop decision (Task J): drop 5 rate_* features (icc_player < 0.20); keep
  `rate_defensive_duels_won` (outfield-only icc_outfield=0.389, keep-tier after GK fix);
  keep `rate_successful_dribbles` (borderline, p99 unremarkable). 48 features retained for rebuild.
- p99 calibration (Task I): expected p99 under t(2.95) = 3.37. Only 6 of 53 features genuinely
  exceed this; 47 are relative leaders among well-behaved residuals.

**Known issues (diagnosed; fixes deferred — full notes in `PROJECT_RUNDOWN_AND_RUNTIME.txt`):**

1. **Runtime is excessive.** K=2 t-model ≈ 8,200 s/chain (~2.3 h wall). Dominant cost is the
   `z_a_raw` uniqueness block (`I×P ≈ 87,450` non-centred params) plus the ν≈3 heavy-tail
   geometry; warmup dominates. Option: marginalize Normal effects (Option 3) — not yet
   prototyped; needs go-ahead.
2. **Player IDs instead of names.** `data/processed/player_map.csv` has only
   `player_index, player_id` — no `player_name`. Plot scripts (12, 13) join then use
   `player_name` *if present*, so they silently fall back to IDs. Fix: rebuild `player_map.csv`
   with names (raw data has `player_name`). Note: `model_objects$metadata$player_name` exists
   and is used by the residual-analysis script (stage 22).
3. **Figure 4 ("PC1 not showing").** The auto-labeler in `12_football_interpretation.R` gives
   **both** factors the same label, so `plot_radar_by_group()` collapses them into one dodged
   series. It also plots raw `Λ`, not the PCA-rotated loadings the caption claims.
4. **Figure 7 grouping.** Season-effect grid facets an arbitrary top-12 by ICC_season; should be
   grouped in a football-meaningful way (reuse `feature_group_lookup()`).
5. **PCA on `Σ_a` vs `ΛΛ'`.** Stage 13 eigendecomposes `ΛΛ'` (common variance only); professor
   suggests `Σ_a = ΛΛ' + Ψ_a`.
6. **Rank K.** No elbow through K=4; Pathfinder screening for K=5–8 is now ready to run
   (`BFA_RANK_METHOD=pathfinder BFA_RANK_MAX=8`); full NUTS to follow on promising ranks.
7. **Low-rank vs diagonal: better metrics added.** Frobenius distance = 4.54, K-fold Δ = +72,881
   (see stage 24). Remaining issue: correlation structure misses key defensive blocs (K=2 too low).
8. **Student-t** ν pinned near its lower bound (≈3). Transform analysis (stage 23) identified 5
   log1p + 18 sqrt candidates. Whether transforms move ν requires a refit (gate: go-ahead for
   Task K rebuild which bundles GK exclusion + feature drops + transforms).
9. **GK players not yet excluded.** Present in the panel as rows; GK feature columns were excluded
   in stage 01, but player rows were not. Fix: `total_gk_saves == 0 | is.na(total_gk_saves)` filter
   in stage 01 before the `eligible` step. Goes into Task K's rebuild bundle.
10. **LOO reliability.** Pareto-k "very bad" for >58% of points. Use K-fold ELPD (already
    computed, Δ = +72,881) as primary metric; PSIS-LOO as secondary with explicit Pareto-k caveat.
11. **Lambda_a convergence — mixing problem escalated.** Raw Lambda_a: min ESS_bulk=19, max
    Rhat=1.185. Lambda_a Lambda_a' (rotation-invariant): diagonal ESS=104, Rhat=1.056;
    off-diagonal ESS=24, Rhat=1.130. 5.5× ESS improvement (rotation-labeling contribution)
    but still below thresholds. GK bimodal structure is the leading hypothesis. Priority: run
    Task K (GK exclusion + rebuild) and check whether ESS improves before proceeding to Task O.
12. **Two-component mixture (Issue 12).** Proposed alternative to Student-t:
    ε ~ π·N(0,σ₁²) + (1-π)·N(0,σ₂²) with σ₁ < σ₂. More interpretable than t for football
    (discrete outlier fraction vs smooth heavy tail). Requires simulation-recovery check before
    fitting real data.

**Next required go-ahead:** Task K (stage 01+ rebuild: GK exclusion + 5 feature drops + 23
feature transforms; then refit stages 10 and 18). Optional: run GK-only ablation fit (stage 18
with only GK exclusion applied, ~2.3h extra) before the full bundle, for clean attribution of
which change moved ν and improved Lambda_a Lambda_a' ESS.

**Chain initialisation — confirmed in place:** `build_pca_init()` (`src/stan_helpers.R:179`) is
wired into stages 10, 18 and 21.

---

## 9. Working agreements

- **Long runs in the background.** Stan fits / pipeline stages are run with
  `run_in_background: true`; report a structured summary on completion. Logs land in `logs/`.
- **Do not resample unless asked** — default to `BFA_REUSE_FIT=true` and reuse existing CSV fits.
- This file and `PROJECT_RUNDOWN_AND_RUNTIME.txt` were produced 2026-06-26 in a documentation-only
  session; the open issues above were diagnosed but intentionally **not** implemented.
