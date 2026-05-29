# Fit low-rank-Sigma_a + diagonal-Sigma_b model to real La Liga data.
# Purpose: discover the actual player factor structure in La Liga performance.
# This is where football interpretation comes from: the Lambda_a loadings
# reveal which combinations of features define latent player archetypes.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/10_fit_real_lowrank_a_diag_b.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

model_objects_path <- file.path(paths$processed, "model_objects.rds")
Y_path             <- file.path(paths$processed, "Y_scaled.csv")
stan_data_path     <- file.path(paths$stan_data, "real_lowrank_a_diag_b_stan_data.rds")

if (!file.exists(model_objects_path)) {
  stop("Missing model_objects.rds. Run stage 02 first.", call. = FALSE)
}
if (!file.exists(Y_path)) {
  stop("Missing Y_scaled.csv. Run stage 02 first.", call. = FALSE)
}

model_objects <- readRDS(model_objects_path)
Y_real        <- as.matrix(readr::read_csv(Y_path, show_col_types = FALSE))

stan_data <- build_lowrank_a_diag_b_stan_data(
  model_objects = model_objects,
  Y             = Y_real,
  rank_a        = simulation$rank_a
)
saveRDS(stan_data, stan_data_path)
message("Stan data saved: ", stan_data_path)

init_fn <- build_pca_init(stan_data)

fit <- fit_cmdstan_model(
  stan_file  = file.path(paths$stan, "additive_lowrank_a_diag_b.stan"),
  stan_data  = stan_data,
  fit_path   = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
  seed       = pipeline$seed + 100L,
  model_id   = "10_real_lowrank_a_diag_b",
  init       = init_fn
)

write_fit_outputs(fit, "10_real_lowrank_a_diag_b")

write_notes(
  c(
    "Real La Liga Data: Low-Rank Sigma_a + Diagonal Sigma_b",
    "",
    paste0("N = ", nrow(Y_real), " player-season rows"),
    paste0("P = ", ncol(Y_real), " scaled features"),
    paste0("I = ", length(model_objects$player_levels), " unique players"),
    paste0("S = ", length(model_objects$season_levels), " seasons"),
    paste0("Q_a = ", simulation$rank_a, " player latent factors"),
    "",
    "Season effects are plain diagonal: Sigma_b = Diag(sigma_b^2).",
    "Player effects are low-rank plus uniqueness: Sigma_a = Lambda_a Lambda_a' + Diag(psi_a^2).",
    "",
    "Loading interpretation in stage 12."
  ),
  file.path(paths$notes, "10_real_lowrank_a_diag_b_notes.txt")
)

message("10_fit_real_lowrank_a_diag_b complete.")
