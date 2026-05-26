# Fit the real-data diagonal additive player-plus-season model.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/03_fit_real_diagonal_additive.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

model_objects_path <- file.path(paths$processed, "model_objects.rds")
if (!file.exists(model_objects_path)) {
  stop("Missing model objects. Run stage 02 first.", call. = FALSE)
}

model_objects <- readRDS(model_objects_path)
stan_data <- build_diagonal_stan_data(model_objects)
saveRDS(stan_data, file.path(paths$stan_data, "real_diagonal_additive_stan_data.rds"))

fit <- fit_cmdstan_model(
  stan_file = file.path(paths$stan, "additive_diagonal.stan"),
  stan_data = stan_data,
  fit_path = file.path(paths$fits, "03_real_diagonal_additive_fit.rds"),
  seed = pipeline$seed + 30L,
  model_id = "03_real_diagonal_additive"
)

write_fit_outputs(fit, "03_real_diagonal_additive")

write_notes(
  c(
    "Real Diagonal Additive Model",
    "",
    "Model: y_ij = mu + a_i + b_j + epsilon_ij.",
    "Player, season, and residual covariance matrices are diagonal in the scaled feature space.",
    "Random effects use non-centered standard-normal parameters and are column-centered in transformed parameters."
  ),
  file.path(paths$notes, "03_real_diagonal_additive_notes.txt")
)

message("03_fit_real_diagonal_additive complete.")
