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
| 22 | `22_residual_icc_screen.R` | Residual analysis + ICC feature screen |
| 23 | `23_transform_recommendations.R` | Transform recommendations per feature |
| 24 | `24_frobenius_kfold_icc.R` | Frobenius dist, K-fold ELPD, ICC shift, Pareto-k |
| 26 | `26_ablation_gk_only_t.R` | GK-only ablation fit (stage 26, K=2 t, N=4586) |
| 27 | `27_sim_marginalized_recovery.R` | Sim recovery for marginalized Normal model (I=200) |
| 28 | `28_fit_k_t_prod.R` | **Production fit:** t-model at chosen K (K=3) |
| B4 | `scripts/B4_postfit_diagnostics.R` | Task L (ν) and Task S (Λ_a mixing) post-fit verdicts |
| B-diag | `scripts/B_block_b_diagnostics_28.R` | Stage 28 diagnostics via Python column-pass |
| B-loo | `scripts/C_block_c_loo_k2_vs_k3.R` | LOO K=2 vs K=3 via generate_quantities |
| B-marg | `scripts/B_block_f_marginalized_i1000.R` | Marginalized sim recovery at I=1000 |
| B-plt | `scripts/B_block_h_loading_plots_k3.R` | K=3 loading plots from stage 28 posterior summary |

Note: stages 22-28 and B-scripts exist but are **not** registered in `run_pipeline.R`'s `stage_map`.

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
  discovery, residual calibration, Lambda_a mixing check, feature decision), Entry 7 (full rebuild:
  GK exclusion, feature drops, transforms, K=3 fit), Entry 8 (Session 9 diagnostics, LOO K=2 vs K=3,
  loading plots, professor summary). **Always read the latest entry when resuming.**
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

## 8. Current state & known issues (as of 2026-08-25)

**State:** **Model near closure — K=3 comparative justification in progress.** Full
rebuild (Tasks K–Q) complete. Production fits: stage 10 (Normal, K=2, P=48, N=4586,
3 chains), stage 18 (Student-t, K=2, P=48, N=4586, 4 chains, ν=4.882), stage 28
(Student-t, **K=3**, P=48, N=4586, 4 chains, ν=4.902), stage 21 (Student-t, K=4, P=48,
N=4586, 4 chains — sampling in progress as of 2026-08-25; launched 2026-08-22).
Minutes-scaled fits (stage 34): variant A (fixed φ=0.5) ν̂=5.667, variant B (estimated
φ) ν̂=5.979, φ̂=0.366. Stage 35 minutes diagnostic: low-minute worst-cell share 8.3%
(14/168 obs) vs 48.9% at stage 28. Verdict: CORRECTED.

**Rank selection (Session 14 extended, 2026-08-25 — per-K comparison):**
K=3 justified from a per-K comparison of K=2 (stage 18), K=3 (stage 28), and K=4
(stage 21, pending) on the three criteria the supervisor prescribed:
(A) variance shares of ΛΛ' PCs — at K=2 last-PC = 30.8% common; at K=3 last-PC = 10.5%;
K=4 last-PC pending (expected near noise floor).
(B) reliable-loader count per PC under fixed-reference PCA — K=2 last-PC 45/48 reliable;
K=3 last-PC 43/48 reliable; K=4 last-PC pending.
(C) residual adequacy under each fit's own t(ν̂) noise model — K=2 mean frac 0.74%,
K=3 mean frac 0.77% (both below the 1% floor, so residuals adequate at both K and NOT a
discriminator between them; K=4 mean frac pending).
CV disclosed honestly as inconclusive across K ∈ {2..6}: stage-32 CSV argmax is K=6,
stage-17 full-range argmax is K=5 (K=3 a local dip); the earlier qmd claim of "both
metrics peak at K=3" was **incorrect** and has been corrected.

**Implemented Session 14 (2026-08-22 initial + 2026-08-25 extension):**

*Initial rewrite (2026-08-22):*
- **§sec-rank rewritten** in `docs/professor/professor_summary.qmd`: three-legged
  qualitative K=3 argument (variance share / reliable-loader count / residual adequacy)
  plus honest CV disclosure paragraph. Callout at line 46 and model-comparison table
  updated accordingly. Open Questions #1 (diagonal-vs-K3 LOO on rebuilt data) marked
  deferred as not blocking closure.
- **§sec-interp** in the qmd: added fixed-reference PCA methods paragraph explicitly
  citing the posterior-mean $\bar{\Lambda}\bar{\Lambda}^\top$ as the anchor
  orientation for CI computation (addresses supervisor's iteration-to-iteration warning).
- **`outputs/notes/32_recommended_k.txt`** deprecated in place — previous single-line
  "3" replaced with a note explaining the CSV argmax is K=6 and K=3 is retained on
  qualitative grounds. Stage-32 CSV and figure kept as archived CV-inconclusiveness
  evidence.

*Extension (2026-08-25) — user noted the initial rewrite argued "K=3 adequate", not "K=3
preferred over neighbouring K". Scope extended to per-K comparison:*
- **`scripts/29b_kcompare_residuals.R`** — K-parameterised residual analysis. Reads a
  fit's posterior summary for σ_e and ν, auto-discovers A/B parameter column indices
  from the CSV header (they differ across K because lambda_a_free size changes), runs
  the same streaming column-pass extractor as stage 29, and appends a row to
  `outputs/tables/29b_kcompare_summary.csv`.
- **`scripts/31b_kcompare_pca.R`** — K-parameterised PCA-with-CI analysis. Same
  mean-rotation approximation as stage 31; appends to
  `outputs/tables/31b_kcompare_{variance_shares,reliable_counts}.csv`.
- **Stage 21 launched** (2026-08-22 16:14 CEST): K=4 NUTS+t fit, `Q_a=4`,
  `max_treedepth=10`, `adapt_delta=0.92`. As of 2026-08-25 chain 1 has entered sampling
  (~50% overall); progress limited by intermittent macOS sleep.
- **§sec-rank rewritten again** — now a comparative argument (K=2 vs K=3, with K=4 row
  auto-appearing when the fit lands). Tables are driven by `read.csv()` on the kcompare
  CSVs, so no qmd edits are needed once stage 21 completes and the K=4 analyses run.
  Final "K=3 verdict" subsection replaces the previous adequacy-framed conclusion.
- **PDF re-render pending**: after K=4 fit completes and 29b/31b run at K=4, one
  `quarto render` call produces the final PDF.

**Implemented Session 13 (2026-08-08):**
- **Stage 34 complete** (`scripts/34_fit_minutes_scaled_t.R`): both variants fit successfully.
  - **34a (fixed φ=0.5):** ν̂=5.667, 95%CI=[5.536, 5.803], ESS_bulk=2,949, Rhat=1.000. ~2.8h wall clock.
  - **34b (estimated φ):** ν̂=5.979, 95%CI=[5.831, 6.126]; φ̂=0.366, 95%CI=[0.360, 0.373],
    ESS_bulk(φ)=4,793, Rhat=1.000. φ̂<0.5: data reject pure Poisson scaling; super-Poisson
    variance (player heterogeneity partially absorbed by factor model) reduces effective exponent.
  - compute_log_lik=0 in both fits; no LOO vs stage 28 available yet.
  - CSV chains: `fits/csv/34a_*/20260804_061259/` (4 × 2.4GB), `fits/csv/34b_*/20260804_090942/`.
  - Posterior summaries: `outputs/tables/34a_*_posterior_summary.csv`, `outputs/tables/34b_*_posterior_summary.csv`.
- **Stage 35 complete** (`scripts/35_minutes_diagnostic_34a.R`): minutes diagnostic on
  stage-34a residuals. Column-pass on stage-34a CSVs (same indices as stage 28, model
  structure identical). Low-minute share (bottom-20th-pct, top-200 worst cells):
  48.9% (stage 28) → 8.3% (14/168 obs, stage 34a). Verdict: CORRECTED.
  Outputs: `outputs/notes/35_minutes_diagnostic_34a.txt`, `outputs/figures/35_*.png`.
- **Documentation updated** — pipeline_journal.qmd Entry 11, CLAUDE.md, professor_summary.qmd,
  professor_summary.pdf rerendered.

**Implemented Session 12 (2026-08-04):**
- **psi_a sigma_floor fix** — `stan/player_season_lowrank_marginalized.stan`: changed
  `vector<lower=0>[P] psi_a;` → `vector<lower=sigma_floor>[P] psi_a;` (matches full models).
  Prevents Σ_a = ΛΛ' + diag(ψ²) from becoming near-singular during optimization. Binary deleted
  and recompiled.
- **Stage 32 — Player-holdout K-selection CV** (`scripts/32_player_holdout_kselect.R`):
  15% player holdout (228 players / 671 obs), stratified by recurrence class. Uses Pathfinder
  on the marginalized Normal model (~1k params; NUTS infeasible due to O(I×P³) gradient cost
  ~150ms/eval × max_depth exploration). Settings: `BFA_PF_PATHS=4 BFA_PF_DRAWS=400 BFA_PF_ITERS=200`.
  Wall clock ~16 min/K → ~2 h for K=1..8. Results **pending** at time of writing.
  - **Phase 1 (base):** constant σ_e, uniform-diagonal Woodbury; governs base K*.
  - **Phase 2 (mv, φ=0.5):** minutes-scaled σ_{n,p} = σ_{e,p}·√(m̄/m_n), rank-1-plus-diagonal
    Woodbury; governs production K* for stage 34. Both metrics from the same Pathfinder draws.
  - Cache: each K's draws matrix saved to `fits/32_marg_k{K}_draws.rds` (not CSV).
  - **Results (K=1..6):** K=3 peaks on base (−34737.2); mv elbow at K=3 (−34211.0). K=4..6
    non-monotonic (Pathfinder instability). K* = 3 confirmed.
  - Pathfinder limitation: all 4 paths use identical init (cmdstanr bug for `pathfinder()`);
    SE=0.0 for K=2..5. K=7..8 not computed (repeated background task kills).
- **Transform language fix** — `professor_summary.qmd`: "Sparse counts" → "Sparse count rates",
  "Moderate Poisson-like counts" → "Moderate-skew count rates". Added paragraph explaining
  features are count rates with Var(rate) ∝ λ/minutes; transforms applied to rates, not raw counts.
- **professor_summary additions** — New §sec-minutes (minutes-scaled variance motivation, formula,
  φ=0.5 justification); updated §sec-kselect-holdout (player-holdout CV description + ELPD table
  chunk); pending analyses updated.

**Implemented Session 10 (2026-07-02/05):**
- **Stage 29 outlier study** (`scripts/29_outlier_study.R`): full standardised-residual analysis on
  stage 28 posterior mean. τ₉₉=3.135 at ν=4.902; 0.77% cells exceed threshold (vs 2.65% expected).
  Top outlier: Juanmi Latasa aerial duels z=16.7. Tables: `29_residual_feature_summary.csv`,
  `29_residual_worst_cells.csv`, `29_residual_skew_triangulation.csv`, `29_player_factor_scores.csv`.
- **Clustered player scatter** (`scripts/29b_clustered_scatter.R`): k-means k=5 on PC1×PC2, with
  PC1-rank-ordered archetype labels (Forwards / Attacking midfielders / Midfielders / Defensive
  midfielders / Centre backs). Verified against Messi, Neymar, Cristiano, Ramos, Xavi.
- **`professor_summary.qmd` rewrite** (all 9 sections): generative model with full LaTeX, data
  pipeline section (GK exclusion, transforms, feature selection), diagnostics, rank selection
  (Pathfinder + PSIS-LOO), factor interpretation (PCA loadings), Student-t section (ν arc),
  outlier study, model comparison table, next steps.
- **Fixes applied during review:** LOO t-vs-Normal table hardcoded (stale CSV); ICC barchart NA
  facet removed; transforms section expanded (log1p/sqrt rationale, alternatives, skewness table);
  clustered scatter archetype labels corrected (PC1 sign inversion); author name fixed to
  Sebastián Serra Peña.

**Implemented Session 9 (2026-07-02):**
- **Block A** — Stage 10 fresh 4-chain refit launched (chain 1 of previous run at 64%, per no-merge rule). In progress.
- **Block B** — Stage 28 full diagnostics via Python column-pass (289 params, 1 batch call/chain). ν=4.902 (ESS=2734, Rhat=1.000). **Task S PASS at K=3:** LLt.1.10 ESS=237>200, Rhat=1.005. Worst Lambda_a Rhat=1.014; 7/144 Lambda_a params with Rhat>1.01. Cached to `outputs/tables/28_diagnostics_cache.csv`.
- **Block C** — LOO K=2 (stage 18) vs K=3 (stage 28) via `generate_quantities(fitted_params = csv_files)`. In progress (resolves 45-min `as_cmdstan_fit()` bottleneck).
- **Block D** — Stage 03 diagonal refit on P=48/N=4586 fresh 4-chain. In progress.
- **Block E** — Stage 20 LOO diagonal vs K=3. Blocked on Block D.
- **Block F** — Marginalized Normal I=1000 sim recovery (5× I=200 baseline); unit test S_i=2 PASS confirmed. In progress.
- **Block G** — P-mismatch audit: 3 Session 8 fixes documented in `outputs/notes/note_p_mismatch_fixes.md`; stale STAGE18_* constants in `src/large_fit_helpers.R` marked.
- **Block H** — K=3 loading plots from stage 28 posterior summary: `28_lambda_a_heatmap.png`, `28_pca_loading_heatmap.png`, `28_factor_group_radar.png`. PCA of ΛΛ': PC1=61.8%, PC2=27.7%, PC3=10.5%.
- **Block H** (professor summary) — `docs/professor/professor_summary.qmd` updated with K=3 model spec, ν=4.90, new stage 28 figures, Pathfinder rank selection table, P=48/N=4586 data section.

**Implemented Session 8 (2026-07-01/02):**
- **Task K (B1-B2):** GK row filter in stage 01 (N: 4944→4586, I: 1650→1529); 5 rate_*
  features dropped from config.R (P: 53→48); 23-feature transform block in stage 02
  (log1p×5, sqrt×18) before Z-scaling.
- **Task L:** ν moved 2.95→4.882 (stage 18) / 4.902 (stage 28). Verdict: transforms
  absorbed heavy tails. t-model still appropriate.
- **Task S (Λ_a mixing):** LLt.1.10 ESS 24→193 (8×), Rhat 1.130→1.005, between-chain SD
  0.0153→0.0032 (5×). Borderline verdict (ESS=193, threshold 200); escalated with flag.
- **B5 rank selection:** Pathfinder K=0-4 calibration → K*=3 (reversal at K=4).
  Full K=0-8 Pathfinder confirmed unreliable at K≥2 (huge inter-run variance).
- **B6 marginalized Normal:** `stan/player_season_lowrank_marginalized.stan` + stage 27
  sim recovery (coverage 78%, math correct).
- **B7 stage 28:** t-model K=3 fit complete (2.4h wall, 4 chains, ν=4.902).
- **B8 post-processing:** stages 11–16, 19 regenerated. Three P-mismatch bugs fixed in
  pipeline scripts (stages 11, 13, 16). Stage 20 deferred (needs new stage 03 on P=48 data).

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

**Implemented Session 6 (2026-06-30):**
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

**Implemented Session 7 (2026-07-01):**
- p99 formula verified: `qt(0.995, df=nu)/sqrt(nu/(nu-2)) = 3.371` is correct for the unit-variance
  standardisation used. Numerical coincidence with `qt(0.99, df=5)=3.365` confirmed as coincidence;
  both thresholds flag the same 6 features. nu=2.9508 hardcoded in script matches cache (2.950752).
- Sign-alignment check on existing stage-18 draws: ZERO column-level sign flips across 4000 draws.
  LT positive-diagonal constraint is working; column-level sign-switching is not the mixing cause.
  Off-diagonal ESS problem is concentrated in one pair: LLt.1.10 (goals×passes, ESS=24).
  Five of seven selected pairs have adequate ESS (228–775). Problem is specific to `Lambda_a[10,1]`
  (per90_passes loading on factor 1), not global.
- GK-only ablation fit (stage 26, `26_real_lowrank_a_diag_b_t_gk_excl`): N=4586, I=1529, ~4684s.
  Results: between-chain SD reduced 1.5× (0.0153→0.0104); LLt.1.10 ESS improved 2.7× (24→64);
  Lambda_a[10,1] ESS=30. GK exclusion is PARTIAL contributor — mixing not resolved. Grand mean
  shifted from −0.452 to −0.614 (outfield-only estimate more interpretable). Full verdict and
  per-chain baseline in `outputs/notes/note_ablation_gk_lambda.md`.
- `write_fit_outputs()` patched: wrapped `fit$summary()` in `tryCatch` to degrade gracefully when
  a custom model_id doesn't match any `summary_variables_for_model` pattern.

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
4. **Do NOT report raw Λ; report PCA of both ΛΛ' and Σ_a.** Víctor: "Λ ni lo reportaría" —
   raw loadings are identified only up to rotation and are uninterpretable. In professor_summary,
   replace raw-Λ heatmap and top-loaders table with CI-filtered PCA loadings from `31_pca_with_ci.R`.
   Also produce Σ_a = ΛΛ' + Ψ_a variant per draw and report eigenvalue shares side by side.
   Session 11 Task 3+4 → `scripts/31_pca_with_ci.R`.
5. **PCA on `Σ_a` vs `ΛΛ'` — covered in Issue 4 (merged).** See above.
6. **Rank K.** ✅ K*=3 confirmed by player-holdout CV (stage 32, K=1..8). Base metric peaks at
   K=3 (ELPD_base=-34737.2, reversal of -26 units at K=4=-34763.4). MV metric shows elbow at K=3
   (-34211.0), K=4 marginal (+1.2 units), K=5 reverses (-34219.4). K=6..8 still running for
   completeness. Production fit at K=3 (stage 28) is confirmed correct. Stage 34 uses K=3.
7. **Low-rank vs diagonal: better metrics added.** Frobenius distance = 4.54, K-fold Δ = +72,881
   (see stage 24). Stage 20 rebuild (K=3 vs diagonal P=48) blocked on Block D (stage 03 refit).
   When complete, output `outputs/tables/20_icc_diagonal_vs_lowrank_rebuilt.csv`.
8. **Student-t ν.** ✅ Resolved — ν moved 2.95→4.902 after transforms; no longer near boundary.
9. **GK players excluded.** ✅ Resolved in Session 8 — 501 rows removed, N=4586, I=1529.
10. **LOO reliability.** Pareto-k "very bad" for >58% of points. Use K-fold ELPD (already
    computed, Δ = +72,881) as primary metric; PSIS-LOO as secondary with explicit Pareto-k caveat.
11. **Lambda_a convergence.** ✅ Resolved (Session 9) — Task S PASS at K=3: LLt.1.10 ESS=237,
    Rhat=1.005. Worst Lambda_a Rhat=1.014; 7/144 entries with Rhat>1.01; acceptable.
12. **Two-component mixture.** Deprioritised — ν≈5 after transforms; Student-t sufficient.
13. **Minutes-in-variance (Issue 13).** ✅ CONFIRMED + FIXED (Sessions 11–13). Stage 30: 48.9%
    of top-200 worst standardised residuals from bottom-20th-percentile of minutes (p<4.4e-18).
    Fix: `ε_{n,p} ~ t(ν, 0, σ_{e,p}·√(m̄/m_n))` (φ=0.5, Poisson-rate theory).
    Stage 34 complete (Session 13): ν̂ rises 4.902→5.667 (34a)/5.979 (34b); φ̂=0.366.
    Stage 35 minutes diagnostic: low-minute share 8.3% (14/168 obs) vs 48.9% at stage 28.
    Verdict: CORRECTED.

**Next required go-ahead (Session 13 — stage 34 + 35 complete):**
- **Stage 34 complete.** ν̂=5.667 (34a, fixed φ=0.5) / 5.979 (34b, estimated φ); φ̂=0.366.
- **Stage 35 complete.** Low-minute worst-cell share: 8.3% (14/168 obs) vs 48.9% at stage 28. Verdict: CORRECTED.
- **Remaining deferred:**
  - Stage 03 diagonal refit → stage 20 LOO → `20_icc_diagonal_vs_lowrank_rebuilt.csv`.
  - Stage 31 (`31_pca_with_ci.R`) — CI-filtered PCA + Σ_a variant for professor_summary.
- Two-component mixture model (Issue 12) deprioritised — ν≈5 after transforms.
- Team-season effects, temporal GP — DEFERRED; close current model first.

**Chain initialisation — confirmed in place:** `build_pca_init()` (`src/stan_helpers.R:179`) is
wired into stages 10, 18, 21, and 28.

---

## 9. Working agreements

- **Long runs in the background.** Stan fits / pipeline stages are run with
  `run_in_background: true`; report a structured summary on completion. Logs land in `logs/`.
- **Do not resample unless asked** — default to `BFA_REUSE_FIT=true` and reuse existing CSV fits.
- This file and `PROJECT_RUNDOWN_AND_RUNTIME.txt` were produced 2026-06-26 in a documentation-only
  session; the open issues above were diagnosed but intentionally **not** implemented.
