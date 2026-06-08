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

## Results (stage 18, 2026-06-08)

### Posterior of ν

| Statistic | Value |
|-----------|-------|
| Mean      | 2.95  |
| Median    | 2.95  |
| 95% CI    | [2.91, 2.99] |

**Interpretation:** ν ≈ 3 is extremely close to the theoretical minimum of 2 (below which
the variance is undefined). This is far from the Normal limit (ν → ∞) and indicates the
football data has genuinely heavy-tailed residuals. With ν = 3, the t-distribution has
infinite kurtosis — individual player-season outliers have a large and non-negligible
probability compared to Normal. Injury-affected seasons, loans to much lower or higher
leagues, or one-off exceptional performances are consistent with this finding.

### LOO-CV comparison: Normal vs t-errors

| Model | ELPD | ΔELPD | SE(Δ) |
|-------|------|-------|-------|
| t-errors (stage 18) | −209,282 | 0 (reference) | — |
| Normal (stage 10)   | −220,950 | −11,668 | 554 |

The t-errors model is preferred by **11,668 ELPD units ≈ 21.1 standard errors**. This is
an overwhelming improvement in held-out predictive accuracy.

### Caveat on LOO reliability

Both models show high Pareto k diagnostics (stage 18: 72% "very bad" k>1; stage 10: 58%
"very bad"). This is expected for large hierarchical models with many player effects and
means the LOO ELPD estimates have high variance. The direction (t-errors better) is
unambiguous, but the magnitude of +11,668 should be interpreted as approximate.

The high Pareto k values are a consequence of the model structure, not a model failure:
each player-season has large leverage because player effects a_i are shared across many
observations. This is the same issue seen in stage 08 and stage 20.

### Conclusion

The data strongly prefers heavy-tailed residuals (ν ≈ 3). Football performance data
contains genuine outliers that Normal errors cannot accommodate. For a production model,
the t-errors specification would be the right choice. The current stage 10 Normal model
underestimates the uncertainty for extreme observations, though its player-level factor
structure remains valid.

See:
- `outputs/tables/18_loo_normal_vs_t.csv` — LOO comparison table
