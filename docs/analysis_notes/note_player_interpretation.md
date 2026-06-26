# Note: interpreting player and season effects

## Why we eigendecompose Λ Λ' rather than use Λ directly

The loading matrix Λ is only identified up to an orthogonal rotation: for any orthogonal Q,
(Λ Q)(Λ Q)' = Λ Λ'. The lower-triangular convention fixes one particular rotation during
sampling, but it is not necessarily the most interpretable one.

The stable object is the player covariance Σ_a = Λ Λ' + Ψ_a. To get interpretable directions
we eigendecompose the low-rank part:

    B_a = Λ Λ' = U D U'     (U: eigenvectors, D: diagonal of eigenvalues, sorted ↓)

The PCA-oriented loadings are:

    L_PCA = U_K D_K^{1/2}

These are the directions of maximum common player variance, ordered by importance. Player scores
on these directions (the "factor scores") are what we plot and cluster.

Stage 13 (`scripts/13_pca_postprocessing.R`) does this for every posterior draw, then summarises
eigenvalues and loadings with posterior means and credible intervals.

## What the two football factors mean (stage 10 real data)

**Factor 1 — "Volume" axis:**
High positive loadings on goals, shots, touches in box, xG. High negative loadings on recoveries,
forward passes, progressive passes. This axis separates pure finishers and box presences
(positive) from ball-winning distributors (negative).

**Factor 2 — "Creation" axis:**
High loadings on assists, key passes, crosses, xG_assist. This captures the assist/creation
dimension orthogonal to raw goal output.

Together the two factors explain the dominant structure of cross-feature correlations in the
player effects, as confirmed by the PCA eigenvalue plot (stage 13): the first two eigenvalues are
clearly larger than the rest.

## Player clustering (stage 14)

Stage 14 (`scripts/14_player_profiles.R`) clusters the 1650 players by their posterior mean
effects â_i (the full 53-dimensional vector). The clustered heatmap in the journal shows groups
that align naturally with playing positions:
- Cluster of forwards: high Factor 1 scores, high goals/shots
- Cluster of central midfielders: near-zero Factor 1, positive creation
- Cluster of defensive midfielders/defenders: negative Factor 1, high recoveries/passes

The PCA scatter (Factor 1 vs Factor 2) provides a player map where position-specific roles
appear as separate clouds.

## Season effects (stage 15)

Stage 15 (`scripts/15_season_profiles.R`) plots b̂_j for each of the 12 seasons (2014/15 –
2025/26). Most features show modest season-level trends; the ICC_season analysis (stage 16)
confirms these are small relative to player-level variation. Exceptions: smart passes and
pressing duels drift noticeably across the 12-year window, reflecting tactical evolution in
La Liga.
