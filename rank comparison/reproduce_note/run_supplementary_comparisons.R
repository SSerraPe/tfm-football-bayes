source("model_comparison.R")

profile = Sys.getenv("PLAYER_MODEL_TASK_PROFILE", "quick")

if (profile == "quick") {
  message("Running quick supplementary-comparison smoke study.")
  scenarios = c(0, 2)
  fit_ranks = 0:3
  n_folds = 2
  sampling_config = default_sampling_config(
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 150,
    iter_sampling = 150,
    refresh = 50
  )
} else if (profile == "thorough") {
  message("Running longer supplementary comparisons for the report.")
  scenarios = c(0, 2)
  fit_ranks = 0:4
  n_folds = 2
  sampling_config = default_sampling_config(
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.99,
    max_treedepth = 14,
    refresh = 250
  )
} else {
  stop("Unknown PLAYER_MODEL_TASK_PROFILE: ", profile)
}

models = compile_player_season_models()
results = list()

for (true_rank in scenarios) {
  message("Supplementary comparisons, truth K = ", true_rank)
  dataset = simulate_player_season_data(
    I = 80,
    J = 6,
    P = 12,
    K_true = true_rank,
    seed = 900 + true_rank
  )

  missing_season = run_bayesian_supplementary_cv(
    dataset = dataset,
    task = "missing_season",
    fit_ranks = fit_ranks,
    n_folds = n_folds,
    fold_seed = 920 + true_rank,
    fit_seed = 3000 + 100 * true_rank,
    compiled_models = models,
    sampling_config = sampling_config
  )
  new_player = run_bayesian_supplementary_cv(
    dataset = dataset,
    task = "new_player",
    fit_ranks = fit_ranks,
    n_folds = n_folds,
    fold_seed = 930 + true_rank,
    fit_seed = 4000 + 100 * true_rank,
    compiled_models = models,
    sampling_config = sampling_config
  )

  results[[as.character(true_rank)]] = list(
    missing_season = missing_season,
    new_player = new_player
  )
}

summary = do.call(
  rbind,
  lapply(names(results), function(true_rank) {
    do.call(
      rbind,
      lapply(names(results[[true_rank]]), function(task) {
        transform(
          results[[true_rank]][[task]]$summary,
          true_rank = as.integer(true_rank),
          task = task
        )
      })
    )
  })
)

selection = do.call(
  rbind,
  lapply(names(results), function(true_rank) {
    do.call(
      rbind,
      lapply(names(results[[true_rank]]), function(task) {
        data.frame(
          true_rank = as.integer(true_rank),
          task = task,
          best_rank = results[[true_rank]][[task]]$best_rank
        )
      })
    )
  })
)

study = list(selection = selection, summary = summary, results = results)
results_dir = file.path("..", "results", "reproduction")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
prefix = paste0("supplementary_comparisons_", profile)
saveRDS(study, file.path(results_dir, paste0(prefix, ".rds")))
write.csv(
  study$selection,
  file.path(results_dir, paste0(prefix, "_selection.csv")),
  row.names = FALSE
)
write.csv(
  study$summary,
  file.path(results_dir, paste0(prefix, "_elpd.csv")),
  row.names = FALSE
)

cat("\nSupplementary-comparison selection:\n")
print(study$selection)
