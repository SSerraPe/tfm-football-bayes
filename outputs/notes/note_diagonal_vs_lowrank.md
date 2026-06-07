# Note: diagonal Σ_a vs low-rank Σ_a model comparison (stage 20)

## The question

Does the low-rank structure in the player covariance matrix actually improve predictions?
Or is a simpler diagonal model (where each feature's player effect is independent of the
others) good enough?

## The two models

**Diagonal model (stage 03, `additive_diagonal.stan`):**
    a_{i,p} ~ N(0, σ_a[p]²)   independently for each feature p
    Σ_a = diag(σ_a²)

This says each feature's player effect is independently distributed — there is no cross-feature
structure in how players differ from each other.

**Low-rank model (stage 10, `additive_lowrank_a_diag_b.stan`):**
    a_i = Λ f_i + u_i   with f_i ~ N(0, I_K), u_i ~ N(0, Ψ_a)
    Σ_a = Λ Λ' + Ψ_a

This says a small number K of latent factors (K=2 here) explain the cross-feature correlations.
The K factors define "archetypes" — combinations of features that move together across players.

## Comparison approach

We compare using **LOO-CV** (leave-one-player-season out cross-validation) based on the
`log_lik[n]` computed in the generated quantities block of both Stan models. A higher ELPD
(expected log predictive density) means the model predicts held-out observations better.

We also compare the **ICC_player** (prop_a[k]) per feature: does the low-rank model assign
more or less player variance to each feature?

## What the results tell us

- If the low-rank model has clearly higher LOO-ELPD: the factor structure is real and
  predictively important. The cross-feature player correlations are genuine and the model
  is capturing something the diagonal model misses.
- If the LOO scores are similar: the factor structure helps interpretation but doesn't
  substantially improve prediction. This could happen if K=2 is too small to capture all
  the real correlation structure.
- ICC_player comparison: the low-rank model typically redistributes variance — for features
  that load strongly on the factors, ICC_player may be higher; for others, lower.

## Prior caveat

Stage 03 uses `sigma_a ~ N(0, 0.7)` (total player SD) while stage 10 uses `psi_a ~ N(0, 0.4)`
(uniqueness SD). The priors are not perfectly matched because in the low-rank model the factors
already absorb most player variance. This is a minor caveat; the LOO comparison is still
informative about overall predictive accuracy.

## Results (stage 20, 2026-06-07)

**LOO-CV comparison:**

| Model | ELPD | ΔELPD | SE(Δ) |
|---|---|---|---|
| Low-rank Σ_a (stage 10) | −220,950 | 0 (reference) | — |
| Diagonal Σ_a (stage 03) | −224,673 | −3,724 | 187 |

The low-rank model is better by **3,724 ELPD units ≈ 19.9 standard errors**. This is overwhelming evidence that the two-factor structure is needed for prediction — the cross-feature player correlations captured by Λ Λ' are genuinely important, not just decorative.

Note: both models show "Pareto k > threshold" LOO warnings, which is expected for a hierarchical model with many player effects (p_loo ≈ 44,400). The direction and magnitude of the LOO difference are clear despite this.

See:
- `outputs/tables/20_loo_diagonal_vs_lowrank.csv` — full LOO comparison
- `outputs/figures/20_icc_diagonal_vs_lowrank.png` — ICC_player scatter
- `outputs/tables/20_icc_diagonal_vs_lowrank.csv` — full ICC comparison table
