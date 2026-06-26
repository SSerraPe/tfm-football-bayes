# Edit only this first section to match the real dataset.
data_file = "REPLACE_WITH_REAL_DATA.csv"
player_column = "REPLACE_WITH_PLAYER_ID_COLUMN"
season_column = "REPLACE_WITH_SEASON_ID_COLUMN"
variable_columns = c(
  "REPLACE_WITH_VARIABLE_1",
  "REPLACE_WITH_VARIABLE_2"
)

# Read the real data and standardize each measurement column.
source("rank_selection.R")
raw_data = read.csv(data_file)
dataset = prepare_player_season_data(
  raw_data,
  player_column,
  season_column,
  variable_columns
)

# Compare the diagonal model with factor models at ranks 1, 2, 3, and 4.
# Five folds are recommended for the real-data analysis.
comparison = run_rank_selection(
  dataset = dataset,
  ranks = 0:4,
  n_folds = 5,
  fold_seed = 11,
  fit_seed = 1000,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.98,
  refresh = 250
)

# Save the two tables needed to interpret the comparison.
dir.create("results", showWarnings = FALSE)
write.csv(comparison$elpd, "results/real_data_elpd.csv", row.names = FALSE)
write.csv(
  comparison$holdouts,
  "results/real_data_holdouts.csv",
  row.names = FALSE
)

cat("\nHeld-out ELPD by fitted rank:\n")
print(comparison$elpd)
cat("\nRank selected by maximum held-out ELPD:\n")
print(comparison$best_rank)
