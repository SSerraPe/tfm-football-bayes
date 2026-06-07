# Note: rank selection via held-out ELPD (stage 17)

## The problem

The factor rank K_a (number of latent player dimensions) is a fixed tuning parameter in Stan —
we cannot sample it directly. We fit a grid of models K ∈ {0, 1, 2, 3, 4} and pick the best one
using predictive accuracy measured on held-out data.

## Why not AIC or BIC?

AIC and BIC penalise model complexity based on the number of parameters, assuming the model is
well-specified. Here the low-rank + diagonal structure is an approximation, so the true number of
"effective parameters" is not well-defined. A held-out predictive criterion measures what we
actually care about: how well does the model predict new observations?

## The cross-validation design

**Unit of cross-validation: (player, variable) pair.**

Instead of holding out entire player-seasons (which would leave some players unobserved in
training), we hold out specific variables for specific players. For example, if pair (player 5,
variable "goals") is assigned to fold 2, then player 5's goals are omitted from training in fold 2,
but all their other variables remain available. This lets the model still estimate player 5's
general effects from the remaining variables.

**Fold construction:**
For n_folds = 5 folds, each variable gets an independent random permutation of the I = 1650
players, then players are assigned to folds in blocks of ~330. This means about 1/5 of (player,
variable) pairs are held out per fold.

**What gets held out:**
For fold m, the held-out set is:

    H_m = { (i, j, k) : player i observed in season j AND (i, k) assigned to fold m }

Each held-out pair (i, k) contributes one held-out value for every season player i was observed.
With 53 variables and 5 folds, roughly 10-11 variables are held out per player per fold, across
all their seasons.

## The ELPD criterion

For each fold m and rank K:

1. Fit model K on training data D \ H_m (all observations not in H_m).
2. For each held-out value y_{ijk}, compute the posterior predictive density averaged over S
   posterior draws:

       p̂(y_{ijk} | D\H_m) = (1/S) Σ_s N(y_{ijk}; a_{i,k}^(s) + b_{j,k}^(s), σ_{e,k}^(s)²)

3. The log score for that value is log p̂(y_{ijk} | D\H_m).

The total ELPD for rank K is:

    ELPD_K = Σ_m Σ_{(i,j,k)∈H_m} log p̂(y_{ijk} | D\H_m)

Higher ELPD is better. We select K* = argmax_K ELPD_K.

## Simulation validation

The professor tested this on 8 simulated datasets (2 per true rank K ∈ {0, 1, 2, 3}). In every
case the maximum ELPD correctly recovered the generating rank. The ELPD differences between
adjacent ranks are clear and consistent — this is a reliable criterion.

## Computational cost

5 folds × 5 ranks = 25 Stan fits. Each fit: 4 chains × (500 warmup + 500 sampling) = 4000
posterior draws. On real data (N=4944, P=53, I=1650) each fit takes several minutes. Total
runtime is of the order of a few hours. Run in the background.
