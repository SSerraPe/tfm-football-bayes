# Player covariance rank comparison

This folder is ready to share. Run R from this directory.

## Analyze the real data

1. Read `methodology_report.pdf`.
2. Edit only the first section of `run_on_real_data.R`.
3. Run `run_on_real_data.R`.

The script compares $K\in\{0,1,2,3,4\}$ using Stan and writes the real-data
tables to `results/`.

## Files

- `run_on_real_data.R`: the only file that needs to be edited.
- `rank_selection.R`: the real-data workflow, including PCA-based starting
  values adapted from the earlier implementation.
- `stan/`: the two Stan models.
- `results/`: saved tables reported in the note and, after fitting, the
  real-data tables.
- `previous_lowrank_plus_diagonal/`: the earlier implementation and its
  documentation, unchanged.
- `reproduce_note/`: optional code that reruns the simulations reported in the
  note.

## What is left out

For each variable, players are divided into folds independently. If
`(player 17, variable 4)` is assigned to a fold, variable 4 is omitted for
player 17 in every observed season. Other variables for player 17 assigned
to the same fold are also omitted in that fit. Variables assigned to other
folds remain available. Thus, several variables for the same player may be
omitted together. With `M` folds, the expected number omitted for one player
in one fitted fold is `p / M`, where `p` is the number of variables.
