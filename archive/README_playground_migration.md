# Migration note: `tfm_playground/` → project structure (2026-09-06)

`tfm_playground/` (the user's original RStudio scratch project) contained the earliest
exploratory work for this thesis: three RMarkdown notebooks, one compiled Stan prototype,
and several cached R objects. This note records what was migrated into the tracked
project structure, what stayed behind and why, and how each superseded piece relates to
the current production pipeline described in `CLAUDE.md`.

## What was migrated

| Original (`tfm_playground/`) | Migrated to | Status |
|---|---|---|
| `similarity_measures.Rmd` (the notebook that actually produced Chapters 2-3's original n=4256/p=23 numbers) | `archive/exploratory_p23_similarity_analysis.R` | Superseded — see below |
| `playground_BMDS.Rmd`'s PCA/FA/CMDS sections | superseded outright by `scripts/02b_descriptive_pca_fa_cmds.R` | Superseded, not separately migrated |
| `playground_BMDS.Rmd`'s Bayesian MDS (BMDS) section | `archive/prototype_bmds_analysis.R` | Prototype, preserved for provenance |
| `player_season_bayes_outputs/player_season_hier_factor.stan` | `stan/prototypes/player_season_hier_factor.stan` | Prototype, preserved for provenance |
| `player_season_bayesian_methodology_full.Rmd` | not migrated | Separate, unrelated prototyping notebook (see below) |

`tfm_playground/` itself was left in place, untouched — it is the user's personal scratch
project and nothing was deleted from it. The large cached artifacts it contains
(`datadistance_obj.Rds`, 138MB; `artifacts/bmds_cosine_subset_N900_p2_*.rds`, 2.9GB;
`.RData`/`.RDataTmp`, ~10GB combined) were **not** copied into the tracked project
structure; `archive/prototype_bmds_analysis.R` loads the BMDS fit directly from its
original location by relative path rather than duplicating it.

## Why Chapters 2-3's numbers changed

`similarity_measures.Rmd` + `playground_BMDS.Rmd` are confirmed (by directly inspecting
the cached `datadistance_obj.Rds`: dim 4256×23, matching feature names) as the real
source of the n=4256/p=23 numbers originally reported in the thesis's Data Source and
Background chapters. That analysis, however, used a smaller, differently-engineered
feature set than the production Bayesian model (`data/processed/model_objects.rds`:
N=4586, P=48). As of this migration, Chapters 2 and 3 (`tfm_latex/data_source.tex`,
`tfm_latex/background.tex`) have been rebased onto the production dataset via
`scripts/02b_descriptive_pca_fa_cmds.R`, so the two chapters and the Bayesian model of
Chapters 4-5 now describe a single, consistent dataset throughout.

`archive/exploratory_p23_similarity_analysis.R` is a cleaned-up, standalone port of the
original p=23 notebook, kept only so the earlier analysis remains independently
reproducible from raw data (it also fixes two bugs found in the original: a `saveRDS`
path that didn't match what was actually loaded, and a lowercase/uppercase filename
mismatch between the two notebooks). It is not used by any current pipeline stage or
thesis chapter.

## The two unrelated prototypes

- **BMDS** (`archive/prototype_bmds_analysis.R`, `playground_BMDS.Rmd`'s BMDS section):
  a completed Bayesian Multidimensional Scaling run (`bayMDS::bmdsMCMC`, cosine
  dissimilarity, a 900-row season-stratified subset since the full sample was
  infeasible for `bayMDS`'s MCMC sampler, p=2). No current thesis chapter reports BMDS
  results — Chapter 3 covers only classical PCA/FA/CMDS, and the production Bayesian
  model (Chapters 4-5) is the additive player+season+residual factor model, not BMDS.
  Kept for provenance in case a future extension of the thesis revisits distance-based
  Bayesian latent-space modelling.
- **Stan prototype** (`stan/prototypes/player_season_hier_factor.stan`, exported by
  `player_season_bayesian_methodology_full.Rmd`): an early, materially simpler draft of
  the same general idea (additive player+season low-rank factors), predating essentially
  every methodological refinement in the production model:
  - Normal likelihood only (no Student-$t$, no degrees-of-freedom parameter).
  - `Lambda_player`/`Lambda_season` are unconstrained, free matrices — **not**
    identified (no LT-PD constraint), unlike the production model's lower-triangular,
    positive-diagonal `Λ`.
  - The season effect is modelled as **low-rank** (`R_season=2`); the production model
    made the opposite, deliberate choice of a **diagonal** season covariance, justified
    by preliminary low-rank fits producing R-hat > 1.70 with only $S=12$ seasons
    (see `tfm_latex/methodology.tex`, §4.4.2).
  - The `transformed parameters` block's `a = a_raw; b = b_raw;` is a no-op — despite
    being introduced as "non-centred" in the surrounding notebook prose, it does not
    actually implement a non-centred parametrisation.
  - No minutes-scaling of the residual variance.
  A compiled binary of this model exists in `tfm_playground/player_season_bayes_outputs/`
  with no corresponding saved fit, suggesting it was compiled but never (successfully,
  or not persistently) sampled.

`player_season_bayesian_methodology_full.Rmd` itself (the notebook that built and
exported the Stan prototype above) was not migrated as a whole: it is a separate,
parallel exploration on a different, 49-feature, all-season feature set, using varimax
(not oblimin) factor rotation and no PCA — it did not produce any of the numbers
currently reported in the thesis and is unrelated to the p=23/P=48 lineage described
above.
