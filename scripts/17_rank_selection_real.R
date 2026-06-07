# Stage 17: Rank selection on real La Liga data via held-out ELPD.
#
# Fits models with K_a in {0, 1, 2, 3, 4} latent player factors using 5-fold
# cross-validation. The fold unit is a (player, variable) pair: when a pair is
# assigned to a fold, that variable is held out for that player in every season.
# The rank with the highest total held-out log predictive density (ELPD) wins.
#
# Output: outputs/rank_selection/real_data_elpd.csv
#         outputs/rank_selection/real_data_holdouts.csv

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/17_rank_selection_real.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- Load the preprocessed real data ----------------------------------------

model_objects_path <- file.path(paths$processed, "model_objects.rds")
Y_path             <- file.path(paths$processed, "Y_scaled.csv")

if (!file.exists(model_objects_path)) stop("Missing model_objects.rds. Run stage 02 first.", call. = FALSE)
if (!file.exists(Y_path))             stop("Missing Y_scaled.csv. Run stage 02 first.", call. = FALSE)

model_objects <- readRDS(model_objects_path)
Y_real        <- as.matrix(read_csv(Y_path, show_col_types = FALSE))

# ---- Build the data frame expected by rank_selection.R ----------------------
# rank_selection.R expects: one row per (player, season) pair, with columns
# for player_id, season_id, and one column per measurement variable.

player_id <- model_objects$player_levels[model_objects$player_index]
season_id <- model_objects$season_levels[model_objects$season_index]

colnames(Y_real) <- model_objects$variable_names

raw_data <- data.frame(
  player_id = player_id,
  season_id = season_id,
  Y_real,
  check.names = FALSE
)

# ---- Source rank_selection.R (changes working dir so Stan paths resolve) -----

rank_comparison_dir <- file.path(model_root, "rank comparison")
if (!dir.exists(rank_comparison_dir)) {
  stop("Cannot find 'rank comparison/' folder at: ", rank_comparison_dir, call. = FALSE)
}

owd <- setwd(rank_comparison_dir)
on.exit(setwd(owd), add = TRUE)

source("rank_selection.R")

# ---- Prepare dataset and run cross-validation --------------------------------

dataset <- prepare_player_season_data(
  raw_data         = raw_data,
  player_column    = "player_id",
  season_column    = "season_id",
  variable_columns = model_objects$variable_names
)

message(
  "Dataset ready: ", dataset$I, " players, ", dataset$J, " seasons, ",
  dataset$P, " variables, ", dataset$N_obs, " observations."
)

comparison <- run_rank_selection(
  dataset          = dataset,
  ranks            = 0:4,
  n_folds          = 5,
  fold_seed        = 11,
  fit_seed         = 1000,
  chains           = 4,
  parallel_chains  = 4,
  iter_warmup      = 500,
  iter_sampling    = 500,
  adapt_delta      = 0.98,
  refresh          = 250
)

# ---- Save outputs (return to original wd first) ------------------------------

setwd(owd)
on.exit(NULL)  # clear the on.exit handler

output_dir <- file.path(model_root, "outputs", "rank_selection")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(comparison$elpd,     file.path(output_dir, "real_data_elpd.csv"),     row.names = FALSE)
write.csv(comparison$holdouts, file.path(output_dir, "real_data_holdouts.csv"), row.names = FALSE)

cat("\nHeld-out ELPD by fitted rank:\n")
print(comparison$elpd)
cat("\nRank selected by maximum held-out ELPD: K* =", comparison$best_rank, "\n")

write_notes(
  c(
    "Stage 17: Rank selection via held-out ELPD",
    "",
    paste0("Dataset: N=", dataset$N_obs, " obs, I=", dataset$I,
           " players, J=", dataset$J, " seasons, P=", dataset$P, " variables"),
    "Folds: 5 (unit = player-variable pair)",
    "Ranks evaluated: K = 0, 1, 2, 3, 4",
    "",
    "Results saved to outputs/rank_selection/",
    paste0("Best rank: K* = ", comparison$best_rank)
  ),
  file.path(paths$notes, "17_rank_selection_notes.txt")
)

message("Stage 17 complete. Best rank: K* = ", comparison$best_rank)
