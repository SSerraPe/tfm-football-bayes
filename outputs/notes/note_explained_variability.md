# Note: explained variability and intraclass correlations (ICC)

## The idea

Each observation y_{n,k} (player-season n, feature k) is decomposed as:

    y_{n,k} = a_{i[n],k} + b_{j[n],k} + ε_{n,k}

The total variance of feature k is:

    Var(y_k) = Var(a_k) + Var(b_k) + Var(ε_k)
             = Σ_a[k,k] + σ_b[k]² + σ_e[k]²

The three shares are:

    prop_a[k]  =  Σ_a[k,k]  / Var(y_k)   ← ICC_player  (how stable across seasons)
    prop_b[k]  =  σ_b[k]²   / Var(y_k)   ← ICC_season  (how much varies league-wide year to year)
    prop_e[k]  =  σ_e[k]²   / Var(y_k)   ← residual unexplained variance

prop_a and prop_b together are the "explained variability" the professor asks for:

    (Σ_a + Σ_b) / (Σ_a + Σ_b + Σ_e)  =  prop_a + prop_b

This is computed in stage 16 and reported as bar charts and scatter plots in the journal.

## Interpretation

Features with high **ICC_player** (prop_a) are stable individual traits — a player's position on
this dimension barely changes from season to season. These are the features the model uses to
characterise player archetypes.

Features with high **ICC_season** (prop_b) reflect league-wide trends — tactical shifts, rule
changes, or external factors that move all players together in a given season.

Features with high **prop_e** (residual) are noisy — either they are genuinely variable across a
player's career, or they are measured with more error.

## Key results (stage 10 real data)

Top features by **ICC_player** (most stable individual traits):
- recoveries:            88.5%
- passes:                86.7%
- touch_in_box:          86.2%
- progressive_passes:    ~85%
- forward_passes:        ~84%

Top features by **ICC_season** (most affected by season-level trends):
- smart_passes:          36.6%
- pressing_duels:        31.4%

Most features have small ICC_season, meaning league-wide yearly shifts are modest relative to
between-player differences. The factor model for player effects (Σ_a low-rank) captures the
strong cross-feature correlations among the high-ICC_player features.
