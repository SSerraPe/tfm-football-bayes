# Fit the diagonal additive model to simulated low-rank data.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/05_fit_sim_diagonal_additive.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

stan_data_path <- file.path(paths$stan_data, "sim_diagonal_additive_stan_data.rds")
if (!file.exists(stan_data_path)) {
  stop("Missing simulated diagonal Stan data. Run stage 04 first.", call. = FALSE)
}

stan_data <- readRDS(stan_data_path)
fit <- fit_cmdstan_model(
  stan_file = file.path(paths$stan, "additive_diagonal.stan"),
  stan_data = stan_data,
  fit_path = file.path(paths$fits, "05_sim_diagonal_additive_fit.rds"),
  seed = pipeline$seed + 50L,
  model_id = "05_sim_diagonal_additive"
)

write_fit_outputs(fit, "05_sim_diagonal_additive")

write_notes(
  c(
    "Simulated Diagonal Additive Model",
    "",
    "This intentionally fits the diagonal additive covariance model to data generated from low-rank plus diagonal player and season effects.",
    "The purpose is misspecification comparison against known simulation covariance diagonals."
  ),
  file.path(paths$notes, "05_sim_diagonal_additive_notes.txt")
)

message("05_fit_sim_diagonal_additive complete.")
