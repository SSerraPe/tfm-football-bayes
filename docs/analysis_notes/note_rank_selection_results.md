# Note: rank selection results (stage 17)

## Selected rank

Among K ∈ {0, 1, 2, 3, 4}, **K = 4** achieves the highest held-out ELPD.

## ELPD table

| Rank K | ELPD (total) | ELPD per held-out value | Δ vs K−1 |
|--------|-------------|------------------------|-----------|
| 0 | −370,699 | −1.415 | — |
| 1 | −328,999 | −1.256 | +41,700 |
| 2 | −297,818 | −1.137 | +31,181 |
| 3 | −278,331 | −1.062 | +19,487 |
| 4 | −262,074 | −1.000 | +16,257 |

Each additional factor adds substantial predictive accuracy. The ELPD per held-out
value improves from −1.415 (diagonal, K=0) to −1.000 (K=4).

## What this means

**K=4 is the best tested rank.** However, the improvements are not levelling off — the
gap from K=3 to K=4 (+16,257) is still large, similar in magnitude to the K=2→K=3 jump
(+19,487). There is no clear "elbow" that would signal the optimal rank within the tested
range. This suggests K ≥ 5 might improve further.

The current main model uses K=2 (stage 10), which was chosen to match the professor's
suggested starting point. The rank selection results show that K=2 captures real structure
(massive improvement over K=0: +72,881) but leaves substantial predictive signal on the
table compared to K=4 (+35,744 from K=2 to K=4).

## Interpretation in football terms

Each factor captures a latent dimension along which players differ. With K=4:
- The first 2 factors likely capture the broad archetypes identified in stage 12
  (attacking output vs defensive/physical metrics).
- Factors 3 and 4 likely capture finer distinctions: position-specific profiles
  (e.g., wide defenders vs central defenders), or style differences within the same
  position role.

The fact that ELPD still improves strongly at K=4 is consistent with the diversity of
player roles in a full La Liga squad (goalkeeper, defenders, midfielders, forwards —
each with multiple sub-types).

## Caveat

The rank comparison uses a separate Stan model family (`player_season_diagonal.stan` /
`player_season_lowrank.stan` from the professor's `rank comparison/` folder) with simpler
priors than the main model. The absolute ELPD values are not directly comparable to the
LOO scores from stages 10 and 20. What matters is the relative ordering across K values,
which is robust.

## Recommendation

A production version of the model should use K=4 (or explore K=5–6). For the current
thesis the K=2 results (stages 10–16) remain valid as a proof of concept: the factor
structure is real and interpretable, and the methodology is sound. Future work should
refit with K ≥ 4.

See:
- `outputs/rank_selection/real_data_elpd.csv` — full ELPD table
- `outputs/rank_selection/real_data_holdouts.csv` — fold structure
- `outputs/notes/note_rank_selection_methodology.md` — how the CV-ELPD procedure works
