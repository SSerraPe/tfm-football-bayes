# Note: the global mean μ and data centering

## What the current model has

The main Stan model (`stan/additive_lowrank_a_diag_b.stan`) declares a `vector[P] mu` parameter
with prior `mu[p] ~ normal(0, 1)` and includes it in the likelihood mean for every observation.
On paper the model is:

    y_n = μ + a_{i[n]} + b_{j[n]} + ε_n

## Why μ is redundant here

Stage 02 (`scripts/02_prepare_model_objects.R`) standardises every feature to mean 0 and standard
deviation 1 before saving `data/processed/Y_scaled.csv`. By the time the data reaches Stan, each
column already has empirical mean ≈ 0.

The player effects a_i and season effects b_j both have priors centred at zero, and the model
enforces within-group centring at every iteration (lines 59–72 of the Stan model). So the only
thing μ absorbs is the (already-removed) overall column mean — its posterior should be very close
to zero for every feature.

The professor's `ideas_and_such/ideas.Rmd` is explicit: *"There is no explicit mean parameter.
The data should therefore be centered in a meaningful way before fitting."* The rank comparison
Stan models (`rank comparison/stan/`) have no μ either.

## Consequence for interpretation

Stage 10 results are not affected: the sampler found μ ≈ 0 and the player/season effects carry all
the interesting signal, exactly as intended. No re-fitting is needed.

## What we do going forward

New Stan models (stages 18 onward) are written without μ, matching the professor's notation:

    y_{ij} = a_i + b_j + ε_{ij}

The data centering step in stage 02 already takes care of the global mean.
