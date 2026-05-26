# Fit low-rank-Sigma_a + diagonal-Sigma_b model to simulated data.
# Purpose: validate that the player factor structure is recoverable when
# Sigma_b is kept diagonal (removing the season factor that caused poor mixing
# in the full low-rank recovery model, stage 06).

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/09_fit_sim_lowrank_a_diag_b.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

model_objects_path <- file.path(paths$processed, "model_objects.rds")
sim_Y_path         <- file.path(paths$processed, "simulated_Y_scaled.csv")
stan_data_path     <- file.path(paths$stan_data, "sim_lowrank_a_diag_b_stan_data.rds")

if (!file.exists(model_objects_path)) {
  stop("Missing model_objects.rds. Run stage 02 first.", call. = FALSE)
}
if (!file.exists(sim_Y_path)) {
  stop("Missing simulated_Y_scaled.csv. Run stage 04 first.", call. = FALSE)
}

model_objects <- readRDS(model_objects_path)
Y_sim         <- as.matrix(readr::read_csv(sim_Y_path, show_col_types = FALSE))

stan_data <- build_lowrank_a_diag_b_stan_data(
  model_objects = model_objects,
  Y             = Y_sim,
  rank_a        = simulation$rank_a
)
saveRDS(stan_data, stan_data_path)
message("Stan data saved: ", stan_data_path)

fit <- fit_cmdstan_model(
  stan_file  = file.path(paths$stan, "additive_lowrank_a_diag_b.stan"),
  stan_data  = stan_data,
  fit_path   = file.path(paths$fits, "09_sim_lowrank_a_diag_b_fit.rds"),
  seed       = pipeline$seed + 90L,
  model_id   = "09_sim_lowrank_a_diag_b"
)

write_fit_outputs(fit, "09_sim_lowrank_a_diag_b")

write_notes(
  c(
    "Simulated Data: Low-Rank Sigma_a + Diagonal Sigma_b",
    "",
    paste0("Player factor rank Q_a = ", simulation$rank_a),
    "Season effects: plain diagonal (sigma_b per feature, no latent factors).",
    "This is the key improvement over stage 06 (full low-rank recovery):",
    "  - Stage 06 had both Sigma_a AND Sigma_b as factor models => poor season mixing (max Rhat 1.739)",
    "  - Stage 09 keeps Sigma_b diagonal => better identifiability with only 12 seasons.",
    "",
    "Recovery diagnostics in stage 11 compare Lambda_a posterior vs simulation truth."
  ),
  file.path(paths$notes, "09_sim_lowrank_a_diag_b_notes.txt")
)

message("09_fit_sim_lowrank_a_diag_b complete.")
