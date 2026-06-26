source("model_comparison.R")

profile = Sys.getenv("PLAYER_MODEL_SIM_PROFILE", "quick")

if (profile == "quick") {
  message("Running quick Bayesian smoke study. Set PLAYER_MODEL_SIM_PROFILE=validation, thorough, or full for longer studies.")
  scenarios = c(0, 2)
  replicates = 1
  fit_ranks = 0:3
  n_folds = 2
  sampling_config = default_sampling_config(
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 150,
    iter_sampling = 150,
    refresh = 50
  )
} else if (profile == "validation") {
  message("Running Bayesian all-ranks validation study.")
  scenarios = 0:3
  replicates = 1
  fit_ranks = 0:4
  n_folds = 2
  sampling_config = default_sampling_config(
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 200,
    iter_sampling = 200,
    refresh = 100
  )
} else if (profile == "thorough") {
  message("Running longer Bayesian replication study for the report.")
  scenarios = 0:3
  replicates = 2
  fit_ranks = 0:4
  n_folds = 2
  sampling_config = default_sampling_config(
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.98,
    refresh = 250
  )
} else if (profile == "full") {
  message("Running full Bayesian replication study.")
  scenarios = c(0, 1, 2, 3)
  replicates = 10
  fit_ranks = 0:4
  n_folds = 5
  sampling_config = default_sampling_config(
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    refresh = 250
  )
} else {
  stop("Unknown PLAYER_MODEL_SIM_PROFILE: ", profile)
}

study = run_bayesian_simulation_study(
  scenarios = scenarios,
  replicates = replicates,
  fit_ranks = fit_ranks,
  I = 80,
  J = 6,
  P = 12,
  n_folds = n_folds,
  seed = 500,
  sampling_config = sampling_config
)

results_dir = file.path("..", "results", "reproduction")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
suffix = profile

saveRDS(
  study,
  file.path(results_dir, paste0("main_simulation_", suffix, ".rds"))
)
write.csv(
  study$runs,
  file.path(results_dir, paste0("main_simulation_runs_", suffix, ".csv")),
  row.names = FALSE
)
write.csv(
  study$cv_summary,
  file.path(results_dir, paste0("main_simulation_elpd_", suffix, ".csv")),
  row.names = FALSE
)
write.csv(
  study$holdouts,
  file.path(results_dir, paste0("main_simulation_holdouts_", suffix, ".csv")),
  row.names = FALSE
)

cat("\nBayesian simulation-study runs:\n")
print(study$runs)
