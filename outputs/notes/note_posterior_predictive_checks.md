# Note: posterior predictive checks (stage 19)

## What are posterior predictive checks?

A posterior predictive check (PPC) asks: if the model is true and we generated new data from
the posterior predictive distribution, would that data look like the observed data? If not,
the model is missing something important.

We run two complementary checks.

## Check 1: cross-variable correlation PPC

**Why:** The key claim of the low-rank model is that the observed cross-feature correlations
among player performance metrics are driven by a small number of latent factors. If this is
true, the model-implied correlation matrix should match the observed sample correlation.

**How:**
- For each of 500 posterior draws, compute the model-implied player covariance:
    Σ_a^(s) = Λ^(s) (Λ^(s))' + diag(ψ_a^(s)²)
- The off-diagonal model-implied covariance between features k and l is Σ_a[k,l], since
  Σ_b and Σ_e are both diagonal (they do not create cross-feature dependence).
- Divide by the total marginal standard deviations to get the model-implied correlation:
    Corr_model(k,l) = Σ_a[k,l] / sqrt(Var(y_k) * Var(y_l))
  where Var(y_k) = Σ_a[k,k] + σ_b[k]^2 + σ_e[k]^2.
- Compare to the observed sample correlation of Y.

**What to look for:** Points close to the diagonal in the scatter plot mean the model is
capturing those feature-pair correlations. Points far from the diagonal indicate covariance
structure that the two-factor model cannot explain (either a third factor is needed, or the
non-Gaussian structure is important).

## Check 2: residual distribution

**Why:** The model assumes residuals ε_{n,p} ~ N(0, σ_e[p]). If player-season performances
have heavy tails or skewness not captured by the player and season effects, the residuals
will not look Normal.

**How:**
- Compute residuals: e_{n,p} = Y[n,p] − â_{i[n],p} − b̂_{j[n],p}
  using posterior mean player effects (from stage 14) and season effects (from stage 15).
- Compare empirical residual histogram to Normal(0, σ_e[p]) for a selection of features:
  the two features with highest ICC_player, the two with highest ICC_season, and two with
  low explained variability.

**What to look for:** If residuals are heavier-tailed than Normal, the Student-t model from
stage 18 should improve fit. If residuals are asymmetric or multimodal, there may be an
unmodelled structure (e.g., position groups).

## Results

See:
- `outputs/figures/19_correlation_ppc_scatter.png` — correlation scatter plot
- `outputs/figures/19_residual_distribution_check.png` — residual histograms
- `outputs/tables/19_correlation_ppc.csv` — full table of observed vs predicted correlations
