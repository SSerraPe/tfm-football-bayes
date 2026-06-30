# Lambda_a Identification Check

## Question

Stage-18 diagnostics (Task F) found raw `Lambda_a` loadings mixing poorly:
min ESS_bulk ≈ 19, max Rhat = 1.185. Does this poor mixing propagate to
`Lambda_a Lambda_a'` (the rotation-invariant quantity underlying all
interpreted thesis outputs: ICC, PCA-rotated loadings)?

## Method

Computed `(Lambda_a Lambda_a')_{pp} = L[p,1]^2 + L[p,2]^2` (diagonal,
P=53 values) and selected off-diagonal pairs from 4000 draws (4 chains × 1000
iterations). ESS_bulk and Rhat computed via `posterior::summarise_draws()`.

## Results

| Quantity | min ESS_bulk | max Rhat |
|---|---|---|
| Raw `Lambda_a` (from Task F) | ~19 | ~1.185 |
| `Lambda_a Lambda_a'` diagonal (P=53) | 104 | 1.056 |
| `Lambda_a Lambda_a'` off-diagonal (7 pairs) | 24 | 1.130 |

ESS improvement factor (diagonal vs raw Lambda_a): 5.5x

**Verdict: escalated (partial rotation-artifact improvement)**

There is a 5× improvement in ESS from raw `Lambda_a` to `Lambda_a Lambda_a'`,
suggesting rotation-labeling is a _contributing_ factor to the poor mixing.
However, both ESS (< 200) and Rhat (> 1.05) remain above acceptable thresholds
for the rotation-invariant product. The posterior geometry itself also has
mixing difficulty.

**Plausible cause:** the bimodal GK/outfield subgroup creates a structural
heterogeneity in the loading space that makes exploration difficult. Removing
GK rows (Task K) is a priority experiment — if it substantially improves
`Lambda_a Lambda_a'` ESS, this supports the GK-geometry hypothesis.

**Recommendation:** do not run Task O at the current K=2 fit without first
running Task K and checking whether GK exclusion resolves the mixing.
If it does not, consider a tighter prior on `lambda_a_free` or Option 3's
marginalized Normal model.

_Written by scripts/25c_lambda_mixing_check.R on 2026-06-30_
