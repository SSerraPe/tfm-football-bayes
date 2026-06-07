# Note: Student-t residuals (stage 18)

## Motivation

The baseline model assumes residuals ε_{n,p} ~ N(0, σ_e[p]). In football data, some
player-seasons may be genuine outliers: injury-affected seasons, loan stints, or simply
extraordinary or catastrophic performances. Normal errors down-weight these too aggressively
and can distort the player effect estimates.

The Student-t distribution with ν degrees of freedom is a heavier-tailed alternative. When
ν → ∞ it converges to Normal; for ν < 10 the tails are substantially heavier and individual
outliers have less influence on the parameter estimates.

## The model change

The Stan model `additive_lowrank_a_diag_b_t.stan` is the same as stage 10 except:
- No global mean μ (data is already Z-scaled — the μ in stage 10 was redundant).
- Likelihood: Y[n,p] ~ student_t(ν, mean_np, σ_e[p]) instead of Normal.
- New parameter: ν ~ Gamma(2, 0.1) (Juárez & Steel 2010; also the Stan reference prior).

The Gamma(2, 0.1) prior puts most mass on ν ∈ [5, 50], allowing moderate to light tails,
while down-weighting ν < 2 (degenerate) and allowing ν → ∞ (Normal).

## How to interpret the results

**Posterior of ν:**
- ν >> 30: t distribution is essentially Normal. The extra parameter is not needed and
  normal errors in stage 10 are fine.
- ν ≈ 10–30: mild heavy tails. Some robustness gain without dramatic changes.
- ν < 10: meaningful heavy tails. Outlier player-seasons are present and the t model
  is materially different from Normal.

**LOO-CV comparison:**
The elpd_diff from `loo_compare()` tells us how much better (or worse) the t model predicts
held-out player-seasons. A positive elpd_diff (normalised by its SE) greater than 2 is
considered clear evidence the t model is preferred.

## Results

See:
- `outputs/tables/18_loo_normal_vs_t.csv` — LOO comparison table
- `outputs/notes/18_real_lowrank_a_diag_b_t_notes.txt` — ν summary after fit completes
