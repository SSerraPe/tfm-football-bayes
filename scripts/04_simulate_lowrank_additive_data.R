# Simulate additive data on the observed row/player/season structure.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/04_simulate_lowrank_additive_data.R"
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
sim <- simulate_lowrank_additive(
  model_objects = model_objects,
  rank_a = simulation$rank_a,
  rank_b = simulation$rank_b,
  seed = simulation$seed
)

sim_model_objects <- model_objects
sim_model_objects$Y <- sim$Y
sim_model_objects$Y_sim <- sim$Y
sim_model_objects$simulation_truth <- sim$truth

diagonal_stan_data <- build_diagonal_stan_data(sim_model_objects, Y = sim$Y)
lowrank_stan_data <- build_lowrank_recovery_stan_data(
  model_objects = sim_model_objects,
  Y = sim$Y,
  rank_a = simulation$rank_a,
  rank_b = simulation$rank_b
)

truth_loadings <- dplyr::bind_rows(
  matrix_to_long_table(sim$truth$Lambda_a, "Lambda_a", model_objects$variable_names, "A"),
  matrix_to_long_table(sim$truth$Lambda_b, "Lambda_b", model_objects$variable_names, "B")
)

truth_variances <- tibble(
  observed_index = seq_along(model_objects$variable_names),
  observed_variable = model_objects$variable_names,
  mu = sim$truth$mu,
  psi_a = sim$truth$Psi_a,
  psi_b = sim$truth$Psi_b,
  sigma_epsilon = sim$truth$sigma_epsilon,
  Sigma_a_diag = sim$truth$Sigma_a_diag,
  Sigma_b_diag = sim$truth$Sigma_b_diag,
  Sigma_epsilon_diag = diag(sim$truth$Sigma_epsilon)
)

simulation_summary <- tibble(
  item = c("seed", "N", "P", "N_players", "N_seasons", "rank_a", "rank_b"),
  value = c(
    simulation$seed,
    nrow(sim$Y),
    ncol(sim$Y),
    length(model_objects$player_levels),
    length(model_objects$season_levels),
    simulation$rank_a,
    simulation$rank_b
  )
)

write_csv(as_tibble(sim$Y, .name_repair = ~ model_objects$variable_names), file.path(paths$processed, "simulated_Y_scaled.csv"))
write_csv(truth_loadings, file.path(paths$tables, "04_simulation_truth_loadings.csv"))
write_csv(truth_variances, file.path(paths$tables, "04_simulation_truth_variance_components.csv"))
write_csv(simulation_summary, file.path(paths$tables, "04_simulation_summary.csv"))
saveRDS(sim_model_objects, file.path(paths$processed, "simulated_lowrank_model_objects.rds"))
saveRDS(sim$truth, file.path(paths$processed, "simulated_lowrank_truth.rds"))
saveRDS(diagonal_stan_data, file.path(paths$stan_data, "sim_diagonal_additive_stan_data.rds"))
saveRDS(lowrank_stan_data, file.path(paths$stan_data, "sim_lowrank_recovery_stan_data.rds"))

message("04_simulate_lowrank_additive_data complete: rank_a=", simulation$rank_a,
        ", rank_b=", simulation$rank_b, ".")
