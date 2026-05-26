# Extensive validation study for the real, simulated-diagonal, and simulated-lowrank fits.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/08_validate_all_models.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

model_objects_path <- file.path(paths$processed, "model_objects.rds")
truth_path <- file.path(paths$processed, "simulated_lowrank_truth.rds")
if (!file.exists(model_objects_path)) {
  stop("Missing model objects. Run stages 00 01 02 first.", call. = FALSE)
}
if (!file.exists(truth_path)) {
  stop("Missing simulated truth. Run stage 04 first.", call. = FALSE)
}

model_objects <- readRDS(model_objects_path)
truth <- readRDS(truth_path)
variable_names <- model_objects$variable_names

model_registry <- tibble(
  model_id = c(
    "03_real_diagonal_additive",
    "05_sim_diagonal_additive",
    "06_sim_lowrank_recovery"
  ),
  label = c(
    "Real data: diagonal additive",
    "Simulated data: diagonal additive misspecification",
    "Simulated data: low-rank recovery"
  ),
  fit_path = file.path(paths$fits, c(
    "03_real_diagonal_additive_fit.rds",
    "05_sim_diagonal_additive_fit.rds",
    "06_sim_lowrank_recovery_fit.rds"
  )),
  summary_path = file.path(paths$tables, c(
    "03_real_diagonal_additive_posterior_summary.csv",
    "05_sim_diagonal_additive_posterior_summary.csv",
    "06_sim_lowrank_recovery_posterior_summary.csv"
  )),
  sampler_path = file.path(paths$diagnostics, c(
    "03_real_diagonal_additive_sampler_diagnostics.csv",
    "05_sim_diagonal_additive_sampler_diagnostics.csv",
    "06_sim_lowrank_recovery_sampler_diagnostics.csv"
  ))
)

artifact_status <- model_registry |>
  mutate(
    csv_manifest_path = file.path(paths$fits, paste0(model_id, "_csv_files.csv")),
    fit_exists = file.exists(fit_path) | file.exists(csv_manifest_path) |
      vapply(model_id, function(id) length(discover_cmdstan_csv_files(id)) > 0L, logical(1)),
    posterior_summary_exists = file.exists(summary_path),
    sampler_diagnostics_exists = file.exists(sampler_path),
    fit_size_mb = ifelse(fit_exists, file.info(fit_path)$size / 1024^2, NA_real_)
  )
write_csv(artifact_status, file.path(paths$diagnostics, "08_validation_artifact_status.csv"))

load_or_create_summary <- function(row) {
  if (file.exists(row$summary_path)) {
    return(read_csv(row$summary_path, show_col_types = FALSE))
  }
  fit <- load_cmdstan_fit(fit_path = row$fit_path, model_id = row$model_id)
  if (is.null(fit)) return(NULL)
  summary_table <- fit$summary(variables = summary_variables_for_model(row$model_id))
  write_csv(summary_table, row$summary_path)
  summary_table
}

summary_list <- setNames(vector("list", nrow(model_registry)), model_registry$model_id)
for (idx in seq_len(nrow(model_registry))) {
  summary_list[[idx]] <- load_or_create_summary(model_registry[idx, ])
}

finite_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

posterior_health <- bind_rows(lapply(names(summary_list), function(model_id) {
  s <- summary_list[[model_id]]
  if (is.null(s)) {
    return(tibble(
      model_id = model_id,
      status = "missing_posterior_summary",
      n_parameters = NA_integer_,
      max_rhat = NA_real_,
      n_rhat_over_1_01 = NA_integer_,
      n_rhat_over_1_05 = NA_integer_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_,
      n_ess_bulk_lt_400 = NA_integer_,
      n_ess_tail_lt_400 = NA_integer_
    ))
  }

  tibble(
    model_id = model_id,
    status = "available",
    n_parameters = nrow(s),
    max_rhat = if ("rhat" %in% names(s)) finite_max(s$rhat) else NA_real_,
    n_rhat_over_1_01 = if ("rhat" %in% names(s)) sum(!is.na(s$rhat) & s$rhat > 1.01) else NA_integer_,
    n_rhat_over_1_05 = if ("rhat" %in% names(s)) sum(!is.na(s$rhat) & s$rhat > 1.05) else NA_integer_,
    min_ess_bulk = if ("ess_bulk" %in% names(s)) finite_min(s$ess_bulk) else NA_real_,
    min_ess_tail = if ("ess_tail" %in% names(s)) finite_min(s$ess_tail) else NA_real_,
    n_ess_bulk_lt_400 = if ("ess_bulk" %in% names(s)) sum(!is.na(s$ess_bulk) & s$ess_bulk < 400) else NA_integer_,
    n_ess_tail_lt_400 = if ("ess_tail" %in% names(s)) sum(!is.na(s$ess_tail) & s$ess_tail < 400) else NA_integer_
  )
}))
write_csv(posterior_health, file.path(paths$diagnostics, "08_posterior_health_summary.csv"))

compute_ebfmi <- function(sampler_df) {
  if (!"energy__" %in% names(sampler_df)) return(NA_real_)
  chain_col <- intersect(c(".chain", "chain"), names(sampler_df))[1]
  if (is.na(chain_col)) {
    energy <- sampler_df$energy__
    return(mean(diff(energy)^2, na.rm = TRUE) / stats::var(energy, na.rm = TRUE))
  }
  sampler_df |>
    group_by(.data[[chain_col]]) |>
    summarise(
      ebfmi = mean(diff(energy__)^2, na.rm = TRUE) / stats::var(energy__, na.rm = TRUE),
      .groups = "drop"
    ) |>
    summarise(min_ebfmi = min(ebfmi, na.rm = TRUE)) |>
    pull(min_ebfmi)
}

sampler_health <- bind_rows(lapply(seq_len(nrow(model_registry)), function(idx) {
  row <- model_registry[idx, ]
  sampler_df <- read_csv_if_exists(row$sampler_path)

  if (is.null(sampler_df)) {
    fit <- load_cmdstan_fit(fit_path = row$fit_path, model_id = row$model_id)
    sampler_df <- tryCatch(
      if (is.null(fit)) NULL else tibble::as_tibble(fit$sampler_diagnostics(format = "df")),
      error = function(e) NULL
    )
    if (!is.null(sampler_df)) write_csv(sampler_df, row$sampler_path)
  }

  if (is.null(sampler_df)) {
    return(tibble(
      model_id = row$model_id,
      status = "missing_sampler_diagnostics",
      n_draws = NA_integer_,
      n_divergent = NA_integer_,
      max_treedepth = NA_real_,
      n_max_treedepth_hits = NA_integer_,
      min_ebfmi = NA_real_,
      mean_accept_stat = NA_real_
    ))
  }

  treedepth_col <- intersect(c("treedepth__", "tree_depth__"), names(sampler_df))[1]
  accept_col <- intersect(c("accept_stat__", "acceptance_stat__"), names(sampler_df))[1]
  tibble(
    model_id = row$model_id,
    status = "available",
    n_draws = nrow(sampler_df),
    n_divergent = if ("divergent__" %in% names(sampler_df)) sum(sampler_df$divergent__ > 0, na.rm = TRUE) else NA_integer_,
    max_treedepth = if (!is.na(treedepth_col)) max(sampler_df[[treedepth_col]], na.rm = TRUE) else NA_real_,
    n_max_treedepth_hits = if (!is.na(treedepth_col)) sum(sampler_df[[treedepth_col]] >= sampler$max_treedepth, na.rm = TRUE) else NA_integer_,
    min_ebfmi = compute_ebfmi(sampler_df),
    mean_accept_stat = if (!is.na(accept_col)) mean(sampler_df[[accept_col]], na.rm = TRUE) else NA_real_
  )
}))
write_csv(sampler_health, file.path(paths$diagnostics, "08_sampler_health_summary.csv"))

recovery_metrics <- function(x, group_cols) {
  x |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      n = n(),
      bias = mean(error, na.rm = TRUE),
      mae = mean(abs(error), na.rm = TRUE),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      median_abs_error = stats::median(abs(error), na.rm = TRUE),
      correlation = if (sum(is.finite(mean) & is.finite(truth_value)) > 1L) {
        stats::cor(mean, truth_value, use = "complete.obs")
      } else {
        NA_real_
      },
      coverage_90 = if (all(c("q5", "q95") %in% names(x))) {
        mean(q5 <= truth_value & q95 >= truth_value, na.rm = TRUE)
      } else {
        NA_real_
      },
      n_truth_outside_90 = if (all(c("q5", "q95") %in% names(x))) {
        sum(!(q5 <= truth_value & q95 >= truth_value), na.rm = TRUE)
      } else {
        NA_integer_
      },
      .groups = "drop"
    )
}

sim_diag_summary <- summary_list[["05_sim_diagonal_additive"]]
if (!is.null(sim_diag_summary)) {
  sim_diag_recovery <- bind_rows(
    extract_vector_summary(sim_diag_summary, "var_a", variable_names) |> mutate(component = "player", truth_value = truth$Sigma_a_diag[index]),
    extract_vector_summary(sim_diag_summary, "var_b", variable_names) |> mutate(component = "season", truth_value = truth$Sigma_b_diag[index]),
    extract_vector_summary(sim_diag_summary, "var_e", variable_names) |> mutate(component = "residual", truth_value = diag(truth$Sigma_epsilon)[index])
  ) |>
    mutate(error = mean - truth_value) |>
    select(component, observed_variable, mean, median, sd, q5, q95, truth_value, error, rhat, ess_bulk, ess_tail)

  write_csv(sim_diag_recovery, file.path(paths$tables, "08_sim_diagonal_recovery_by_feature.csv"))
  write_csv(
    recovery_metrics(sim_diag_recovery, "component"),
    file.path(paths$tables, "08_sim_diagonal_recovery_metrics.csv")
  )
}

sim_lr_summary <- summary_list[["06_sim_lowrank_recovery"]]
if (!is.null(sim_lr_summary)) {
  lambda_a_summary <- extract_matrix_summary(sim_lr_summary, "Lambda_a", variable_names) |>
    mutate(truth_value = truth$Lambda_a[cbind(observed_index, factor_index)])
  lambda_b_summary <- extract_matrix_summary(sim_lr_summary, "Lambda_b", variable_names) |>
    mutate(truth_value = truth$Lambda_b[cbind(observed_index, factor_index)])
  loading_recovery <- bind_rows(lambda_a_summary, lambda_b_summary) |>
    mutate(error = mean - truth_value) |>
    select(parameter, observed_variable, factor_index, mean, median, sd, q5, q95, truth_value, error, rhat, ess_bulk, ess_tail)

  variance_recovery <- bind_rows(
    extract_vector_summary(sim_lr_summary, "psi_a", variable_names) |>
      mutate(parameter = "psi_a_sd", truth_value = sqrt(truth$Psi_a[index])),
    extract_vector_summary(sim_lr_summary, "psi_b", variable_names) |>
      mutate(parameter = "psi_b_sd", truth_value = sqrt(truth$Psi_b[index])),
    extract_vector_summary(sim_lr_summary, "sigma_e", variable_names) |>
      mutate(parameter = "sigma_epsilon_sd", truth_value = truth$sigma_epsilon[index]),
    extract_vector_summary(sim_lr_summary, "Sigma_a_diag", variable_names) |>
      mutate(parameter = "Sigma_a_diag", truth_value = truth$Sigma_a_diag[index]),
    extract_vector_summary(sim_lr_summary, "Sigma_b_diag", variable_names) |>
      mutate(parameter = "Sigma_b_diag", truth_value = truth$Sigma_b_diag[index])
  ) |>
    mutate(error = mean - truth_value) |>
    select(parameter, observed_variable, mean, median, sd, q5, q95, truth_value, error, rhat, ess_bulk, ess_tail)

  write_csv(loading_recovery, file.path(paths$tables, "08_sim_lowrank_loading_recovery_by_parameter.csv"))
  write_csv(variance_recovery, file.path(paths$tables, "08_sim_lowrank_variance_recovery_by_feature.csv"))
  write_csv(
    recovery_metrics(loading_recovery, "parameter"),
    file.path(paths$tables, "08_sim_lowrank_loading_recovery_metrics.csv")
  )
  write_csv(
    recovery_metrics(variance_recovery, "parameter"),
    file.path(paths$tables, "08_sim_lowrank_variance_recovery_metrics.csv")
  )
}

compute_loo_summary <- function(row) {
  run_loo <- tolower(Sys.getenv("BFA_RUN_LOO", unset = "false")) %in% c("true", "1", "yes", "y")
  if (!run_loo) {
    return(tibble(model_id = row$model_id, status = "skipped_set_BFA_RUN_LOO_true"))
  }

  fit <- load_cmdstan_fit(fit_path = row$fit_path, model_id = row$model_id)
  if (is.null(fit)) {
    return(tibble(model_id = row$model_id, status = "missing_fit"))
  }
  if (!requireNamespace("loo", quietly = TRUE) || !requireNamespace("posterior", quietly = TRUE)) {
    return(tibble(model_id = row$model_id, status = "missing_loo_or_posterior_package"))
  }

  log_lik <- tryCatch(
    posterior::as_draws_matrix(fit$draws("log_lik")),
    error = function(e) NULL
  )
  if (is.null(log_lik)) {
    return(tibble(model_id = row$model_id, status = "missing_log_lik_draws"))
  }

  loo_fit <- loo::loo(log_lik)
  estimates <- as.data.frame(loo_fit$estimates)
  tibble(
    model_id = row$model_id,
    status = "available",
    elpd_loo = estimates["elpd_loo", "Estimate"],
    se_elpd_loo = estimates["elpd_loo", "SE"],
    p_loo = estimates["p_loo", "Estimate"],
    looic = estimates["looic", "Estimate"],
    se_looic = estimates["looic", "SE"],
    pareto_k_gt_0_7 = sum(loo::pareto_k_values(loo_fit) > 0.7, na.rm = TRUE),
    pareto_k_gt_1 = sum(loo::pareto_k_values(loo_fit) > 1, na.rm = TRUE),
    max_pareto_k = max(loo::pareto_k_values(loo_fit), na.rm = TRUE)
  )
}

loo_summary <- bind_rows(lapply(seq_len(nrow(model_registry)), function(idx) {
  compute_loo_summary(model_registry[idx, ])
}))
write_csv(loo_summary, file.path(paths$tables, "08_loo_summary.csv"))

if ("ggplot2" %in% rownames(installed.packages())) {
  suppressPackageStartupMessages(library(ggplot2))

  if (nrow(filter(posterior_health, status == "available")) > 0L) {
    p_health <- posterior_health |>
      filter(status == "available") |>
      select(model_id, max_rhat, min_ess_bulk, min_ess_tail) |>
      pivot_longer(-model_id, names_to = "metric", values_to = "value") |>
      ggplot(aes(model_id, value, fill = metric)) +
      geom_col(position = "dodge") +
      coord_flip() +
      labs(x = NULL, y = "Value", fill = "Metric") +
      theme_minimal(base_size = 10)
    ggsave(file.path(paths$figures, "08_posterior_health_summary.png"), p_health, width = 8, height = 4.5, dpi = 150)
  }

  if (file.exists(file.path(paths$tables, "08_sim_diagonal_recovery_by_feature.csv"))) {
    sim_diag_recovery <- read_csv(file.path(paths$tables, "08_sim_diagonal_recovery_by_feature.csv"), show_col_types = FALSE)
    p_diag <- ggplot(sim_diag_recovery, aes(truth_value, mean, color = component)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.8) +
      labs(x = "Truth", y = "Posterior mean", color = "Component") +
      theme_minimal(base_size = 10)
    ggsave(file.path(paths$figures, "08_sim_diagonal_recovery.png"), p_diag, width = 7, height = 5, dpi = 150)
  }

  if (file.exists(file.path(paths$tables, "08_sim_lowrank_loading_recovery_by_parameter.csv"))) {
    loading_recovery <- read_csv(file.path(paths$tables, "08_sim_lowrank_loading_recovery_by_parameter.csv"), show_col_types = FALSE)
    p_load <- ggplot(loading_recovery, aes(truth_value, mean, color = parameter)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.75) +
      labs(x = "Truth", y = "Posterior mean", color = "Matrix") +
      theme_minimal(base_size = 10)
    ggsave(file.path(paths$figures, "08_sim_lowrank_loading_recovery.png"), p_load, width = 7, height = 5, dpi = 150)
  }
}

available_fits <- sum(artifact_status$fit_exists)
available_summaries <- sum(artifact_status$posterior_summary_exists)
report_lines <- c(
  "# Model Validation Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Artifact Status",
  "",
  paste0("- Fit files present: ", available_fits, " / ", nrow(artifact_status)),
  paste0("- Posterior summaries present: ", available_summaries, " / ", nrow(artifact_status)),
  "",
  "## Posterior Health",
  "",
  paste(capture.output(print(posterior_health, n = Inf)), collapse = "\n"),
  "",
  "## Sampler Health",
  "",
  paste(capture.output(print(sampler_health, n = Inf)), collapse = "\n"),
  "",
  "## Simulation Recovery Outputs",
  "",
  "- Diagonal misspecification recovery tables are written when `05_sim_diagonal_additive` posterior summaries exist.",
  "- Low-rank loading and variance recovery tables are written when `06_sim_lowrank_recovery` posterior summaries exist.",
  "- LOO summaries are written when fit objects and the `loo` package are available.",
  "",
  "Primary CSV outputs:",
  "",
  "- `outputs/diagnostics/08_posterior_health_summary.csv`",
  "- `outputs/diagnostics/08_sampler_health_summary.csv`",
  "- `outputs/tables/08_sim_diagonal_recovery_metrics.csv`",
  "- `outputs/tables/08_sim_lowrank_loading_recovery_metrics.csv`",
  "- `outputs/tables/08_sim_lowrank_variance_recovery_metrics.csv`",
  "- `outputs/tables/08_loo_summary.csv`"
)

write_notes(report_lines, file.path(paths$notes, "08_model_validation_report.md"))

message("08_validate_all_models complete: fit files present ", available_fits, " / ", nrow(artifact_status), ".")
