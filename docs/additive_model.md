# Additive Bayesian Model

For observed player-season row `n`, player `i[n]`, season `j[n]`, and `P` scaled features:

```text
y_n = mu + a_i[n] + b_j[n] + epsilon_n
```

The real-data model uses diagonal covariance components:

```text
a_i ~ N(0, diag(sigma_a^2))
b_j ~ N(0, diag(sigma_b^2))
epsilon_n ~ N(0, diag(sigma_e^2))
```

Player and season effects are non-centered and column-centered so the intercept `mu` remains identifiable.

The simulation keeps the same observed rows, player indices, season indices, and selected feature names, but generates player and season effects from low-rank plus diagonal covariance:

```text
Sigma_a = Lambda_a Lambda_a' + diag(Psi_a)
Sigma_b = Lambda_b Lambda_b' + diag(Psi_b)
```

Defaults are `rank_a = 2` and `rank_b = 1`, with `rank_a > rank_b`. The low-rank recovery model uses pure-anchor loading constraints: the first `Q` features are positive own-factor anchors and have zero off-factor loadings.

Recovery targets are `Lambda_a`, `Lambda_b`, `Psi_a`, `Psi_b`, `Sigma_a` diagonal, `Sigma_b` diagonal, residual scales, and feature means.
