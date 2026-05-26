# Shared CmdStan and posterior-summary helpers.

build_diagonal_stan_data <- function(model_objects, Y = model_objects$Y) {
  list(
    N = nrow(Y),
    P = ncol(Y),
    I = length(model_objects$player_levels),
    S = length(model_objects$season_levels),
    Y = unname(Y),
    player_index = as.integer(model_objects$player_index),
    season_index = as.integer(model_objects$season_index),
    sigma_floor = pipeline$sigma_floor
  )
}

build_lowrank_recovery_stan_data <- function(model_objects, Y, rank_a, rank_b) {
  P <- ncol(Y)
  list(
    N = nrow(Y),
    P = P,
    I = length(model_objects$player_levels),
    S = length(model_objects$season_levels),
    Q_a = rank_a,
    Q_b = rank_b,
    N_lambda_a_free = anchor_free_count(P, rank_a),
    N_lambda_b_free = anchor_free_count(P, rank_b),
    Y = unname(Y),
    player_index = as.integer(model_objects$player_index),
    season_index = as.integer(model_objects$season_index),
    lambda_a_anchor = default_anchor_index(P, rank_a),
    lambda_b_anchor = default_anchor_index(P, rank_b),
    sigma_floor = pipeline$sigma_floor
  )
}

build_lowrank_a_diag_b_stan_data <- function(model_objects, Y, rank_a) {
  P <- ncol(Y)
  list(
    N = nrow(Y),
    P = P,
    I = length(model_objects$player_levels),
    S = length(model_objects$season_levels),
    Q_a = rank_a,
    N_lambda_a_free = anchor_free_count(P, rank_a),
    Y = unname(Y),
    player_index = as.integer(model_objects$player_index),
    season_index = as.integer(model_objects$season_index),
    lambda_a_anchor = default_anchor_index(P, rank_a),
    sigma_floor = pipeline$sigma_floor
  )
}

fit_csv_manifest_path <- function(model_id) {
  file.path(paths$fits, paste0(model_id, "_csv_files.csv"))
}

discover_cmdstan_csv_files <- function(model_id) {
  csv_root <- file.path(paths$fits, "csv", model_id)
  if (!dir.exists(csv_root)) return(character(0))

  run_dirs <- list.dirs(csv_root, recursive = FALSE, full.names = TRUE)
  if (length(run_dirs) == 0L) return(character(0))
  run_dirs <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]

  for (run_dir in run_dirs) {
    csv_files <- sort(list.files(run_dir, pattern = "[.]csv$", full.names = TRUE))
    if (length(csv_files) > 0L && all(file.exists(csv_files))) {
      return(csv_files)
    }
  }

  character(0)
}

persist_fit_pointer <- function(model_id, stan_file, csv_files, fit_path) {
  readr::write_csv(
    tibble::tibble(model_id = model_id, csv_file = csv_files),
    fit_csv_manifest_path(model_id)
  )
  saveRDS(
    list(
      model_id = model_id,
      stan_file = stan_file,
      csv_files = csv_files,
      created_at = Sys.time(),
      note = "Lightweight fit pointer; reconstruct with cmdstanr::as_cmdstan_fit(csv_files)."
    ),
    fit_path
  )
}

load_cmdstan_fit <- function(fit_path = NULL, model_id = NULL) {
  check_packages(c("cmdstanr"))

  if (!is.null(fit_path) && file.exists(fit_path)) {
    fit_object <- readRDS(fit_path)
    if (inherits(fit_object, "CmdStanMCMC")) {
      return(fit_object)
    }
    if (is.list(fit_object) && !is.null(fit_object$csv_files)) {
      return(cmdstanr::as_cmdstan_fit(fit_object$csv_files))
    }
  }

  if (!is.null(model_id)) {
    manifest_path <- fit_csv_manifest_path(model_id)
    if (file.exists(manifest_path)) {
      manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
      csv_files <- manifest$csv_file[file.exists(manifest$csv_file)]
      if (length(csv_files) > 0L) {
        return(cmdstanr::as_cmdstan_fit(csv_files))
      }
    }

    csv_files <- discover_cmdstan_csv_files(model_id)
    if (length(csv_files) > 0L) {
      return(cmdstanr::as_cmdstan_fit(csv_files))
    }
  }

  NULL
}

summary_variables_for_model <- function(model_id) {
  if (grepl("lowrank_recovery", model_id)) {
    return(c(
      "lp__", "mu", "lambda_a_diag", "lambda_a_free", "lambda_b_diag",
      "lambda_b_free", "Lambda_a", "Lambda_b", "psi_a", "psi_b", "sigma_e",
      "Sigma_a_diag", "Sigma_b_diag", "Sigma_epsilon_diag"
    ))
  }

  if (grepl("lowrank_a_diag_b", model_id)) {
    return(c(
      "lp__", "mu", "lambda_a_diag", "lambda_a_free", "Lambda_a",
      "psi_a", "sigma_b", "sigma_e",
      "Sigma_a_diag", "var_b", "var_e", "prop_a", "prop_b", "prop_e"
    ))
  }

  c(
    "lp__", "mu", "sigma_a", "sigma_b", "sigma_e",
    "var_a", "var_b", "var_e", "prop_a", "prop_b", "prop_e"
  )
}

fit_cmdstan_model <- function(stan_file, stan_data, fit_path, seed, model_id) {
  if (sampler$reuse_fit) {
    reusable_fit <- load_cmdstan_fit(fit_path = fit_path, model_id = model_id)
    if (!is.null(reusable_fit)) {
      message("Reusing existing CmdStan CSV fit for ", model_id)
      persist_fit_pointer(model_id, stan_file, reusable_fit$output_files(), fit_path)
      return(reusable_fit)
    }
  }

  if (!sampler$run_stan) {
    message("Stan sampling skipped for ", model_id,
            ". Set BFA_RUN_STAN=true to sample.")
    return(NULL)
  }

  check_packages(c("cmdstanr"))
  model <- cmdstanr::cmdstan_model(stan_file)
  output_dir <- file.path(paths$fits, "csv", model_id, format(Sys.time(), "%Y%m%d_%H%M%S"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- model$sample(
    data = stan_data,
    seed = seed,
    chains = sampler$chains,
    parallel_chains = sampler$parallel_chains,
    iter_warmup = sampler$iter_warmup,
    iter_sampling = sampler$iter_sampling,
    adapt_delta = sampler$adapt_delta,
    max_treedepth = sampler$max_treedepth,
    refresh = sampler$refresh,
    output_dir = output_dir,
    output_basename = model_id
  )

  csv_files <- fit$output_files()
  persist_fit_pointer(model_id, stan_file, csv_files, fit_path)
  fit
}

write_fit_outputs <- function(fit, model_id) {
  if (is.null(fit)) {
    readr::write_csv(
      tibble::tibble(model_id = model_id, status = "stan_skipped"),
      file.path(paths$diagnostics, paste0(model_id, "_diagnostic_summary.csv"))
    )
    return(invisible(NULL))
  }

  posterior_summary <- fit$summary(variables = summary_variables_for_model(model_id))
  readr::write_csv(
    posterior_summary,
    file.path(paths$tables, paste0(model_id, "_posterior_summary.csv"))
  )

  sampler_diagnostics <- tryCatch(
    fit$sampler_diagnostics(format = "df"),
    error = function(e) NULL
  )
  if (!is.null(sampler_diagnostics)) {
    readr::write_csv(
      tibble::as_tibble(sampler_diagnostics),
      file.path(paths$diagnostics, paste0(model_id, "_sampler_diagnostics.csv"))
    )
  }

  diagnostic_summary <- posterior_summary |>
    dplyr::summarise(
      model_id = model_id,
      max_rhat = max(rhat, na.rm = TRUE),
      min_ess_bulk = min(ess_bulk, na.rm = TRUE),
      min_ess_tail = min(ess_tail, na.rm = TRUE),
      n_rhat_over_1_01 = sum(!is.na(rhat) & rhat > 1.01)
    )
  readr::write_csv(
    diagnostic_summary,
    file.path(paths$diagnostics, paste0(model_id, "_diagnostic_summary.csv"))
  )

  invisible(posterior_summary)
}
