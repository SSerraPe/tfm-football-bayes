# Fit the low-rank recovery model to simulated data.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/06_fit_sim_lowrank_recovery.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

stan_data_path <- file.path(paths$stan_data, "sim_lowrank_recovery_stan_data.rds")
if (!file.exists(stan_data_path)) {
  stop("Missing simulated low-rank Stan data. Run stage 04 first.", call. = FALSE)
}

stan_data <- readRDS(stan_data_path)
fit <- fit_cmdstan_model(
  stan_file = file.path(paths$stan, "additive_lowrank_recovery.stan"),
  stan_data = stan_data,
  fit_path = file.path(paths$fits, "06_sim_lowrank_recovery_fit.rds"),
  seed = pipeline$seed + 60L,
  model_id = "06_sim_lowrank_recovery"
)

write_fit_outputs(fit, "06_sim_lowrank_recovery")

write_notes(
  c(
    "Simulated Low-Rank Recovery Model",
    "",
    "This model matches the low-rank plus diagonal simulation design.",
    "Loadings use pure-anchor constraints: the first Q variables are positive own-factor anchors and have zero off-factor loadings.",
    "Recovery diagnostics compare posterior summaries against saved Lambda, Psi, and covariance-diagonal truth."
  ),
  file.path(paths$notes, "06_sim_lowrank_recovery_notes.txt")
)

message("06_fit_sim_lowrank_recovery complete.")
