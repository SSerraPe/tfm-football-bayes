# Produce compact diagnostics, recovery tables, and figures from saved fits.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/07_diagnostics_and_recovery_plots.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

model_objects_path <- file.path(paths$processed, "model_objects.rds")
truth_path <- file.path(paths$processed, "simulated_lowrank_truth.rds")
if (!file.exists(model_objects_path)) {
  stop("Missing model objects. Run stage 02 first.", call. = FALSE)
}
if (!file.exists(truth_path)) {
  stop("Missing simulated truth. Run stage 04 first.", call. = FALSE)
}

model_objects <- readRDS(model_objects_path)
truth <- readRDS(truth_path)
variable_names <- model_objects$variable_names

real_summary <- read_csv_if_exists(file.path(paths$tables, "03_real_diagonal_additive_posterior_summary.csv"))
sim_diag_summary <- read_csv_if_exists(file.path(paths$tables, "05_sim_diagonal_additive_posterior_summary.csv"))
sim_lr_summary <- read_csv_if_exists(file.path(paths$tables, "06_sim_lowrank_recovery_posterior_summary.csv"))

diagnostic_availability <- tibble(
  artifact = c("real_diagonal_summary", "sim_diagonal_summary", "sim_lowrank_summary"),
  available = c(!is.null(real_summary), !is.null(sim_diag_summary), !is.null(sim_lr_summary))
)
write_csv(diagnostic_availability, file.path(paths$diagnostics, "07_diagnostic_artifact_availability.csv"))

empirical_real <- empirical_variance_decomposition(model_objects)
write_csv(empirical_real, file.path(paths$tables, "07_real_empirical_variance_decomposition.csv"))

if (!is.null(real_summary)) {
  real_var <- dplyr::bind_rows(
    extract_vector_summary(real_summary, "var_a", variable_names) |> mutate(component = "player"),
    extract_vector_summary(real_summary, "var_b", variable_names) |> mutate(component = "season"),
    extract_vector_summary(real_summary, "var_e", variable_names) |> mutate(component = "residual")
  ) |>
    select(component, observed_variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

  real_var_wide <- real_var |>
    select(component, observed_variable, mean) |>
    pivot_wider(names_from = component, values_from = mean, names_prefix = "posterior_var_") |>
    mutate(posterior_var_total = posterior_var_player + posterior_var_season + posterior_var_residual)

  write_csv(real_var, file.path(paths$tables, "07_real_posterior_variance_components.csv"))
  write_csv(real_var_wide, file.path(paths$tables, "07_real_variance_decomposition.csv"))
} else {
  write_csv(empirical_real, file.path(paths$tables, "07_real_variance_decomposition.csv"))
}

truth_variance <- tibble(
  observed_index = seq_along(variable_names),
  observed_variable = variable_names,
  true_var_a = truth$Sigma_a_diag,
  true_var_b = truth$Sigma_b_diag,
  true_var_e = diag(truth$Sigma_epsilon),
  true_psi_a = truth$Psi_a,
  true_psi_b = truth$Psi_b,
  true_sigma_epsilon = truth$sigma_epsilon
)

if (!is.null(sim_diag_summary)) {
  sim_diag_recovery <- dplyr::bind_rows(
    extract_vector_summary(sim_diag_summary, "var_a", variable_names) |> mutate(component = "player", truth_value = truth$Sigma_a_diag[index]),
    extract_vector_summary(sim_diag_summary, "var_b", variable_names) |> mutate(component = "season", truth_value = truth$Sigma_b_diag[index]),
    extract_vector_summary(sim_diag_summary, "var_e", variable_names) |> mutate(component = "residual", truth_value = diag(truth$Sigma_epsilon)[index])
  ) |>
    mutate(error = mean - truth_value) |>
    select(component, observed_variable, mean, median, sd, q5, q95, truth_value, error, rhat, ess_bulk, ess_tail)
  write_csv(sim_diag_recovery, file.path(paths$tables, "07_sim_diagonal_variance_vs_truth.csv"))
}

if (!is.null(sim_lr_summary)) {
  lambda_a_summary <- extract_matrix_summary(sim_lr_summary, "Lambda_a", variable_names) |>
    mutate(truth_value = truth$Lambda_a[cbind(observed_index, factor_index)])
  lambda_b_summary <- extract_matrix_summary(sim_lr_summary, "Lambda_b", variable_names) |>
    mutate(truth_value = truth$Lambda_b[cbind(observed_index, factor_index)])
  loading_recovery <- bind_rows(lambda_a_summary, lambda_b_summary) |>
    mutate(error = mean - truth_value) |>
    select(parameter, observed_variable, factor_index, mean, median, sd, q5, q95, truth_value, error, rhat, ess_bulk, ess_tail)

  psi_a_summary <- extract_vector_summary(sim_lr_summary, "psi_a", variable_names) |>
    mutate(parameter = "psi_a_sd", truth_value = sqrt(truth$Psi_a[index]))
  psi_b_summary <- extract_vector_summary(sim_lr_summary, "psi_b", variable_names) |>
    mutate(parameter = "psi_b_sd", truth_value = sqrt(truth$Psi_b[index]))
  sigma_e_summary <- extract_vector_summary(sim_lr_summary, "sigma_e", variable_names) |>
    mutate(parameter = "sigma_epsilon_sd", truth_value = truth$sigma_epsilon[index])
  Sigma_a_summary <- extract_vector_summary(sim_lr_summary, "Sigma_a_diag", variable_names) |>
    mutate(parameter = "Sigma_a_diag", truth_value = truth$Sigma_a_diag[index])
  Sigma_b_summary <- extract_vector_summary(sim_lr_summary, "Sigma_b_diag", variable_names) |>
    mutate(parameter = "Sigma_b_diag", truth_value = truth$Sigma_b_diag[index])

  variance_recovery <- bind_rows(
    psi_a_summary,
    psi_b_summary,
    sigma_e_summary,
    Sigma_a_summary,
    Sigma_b_summary
  ) |>
    mutate(error = mean - truth_value) |>
    select(parameter, observed_variable, mean, median, sd, q5, q95, truth_value, error, rhat, ess_bulk, ess_tail)

  write_csv(loading_recovery, file.path(paths$tables, "07_sim_lowrank_recovery_loadings.csv"))
  write_csv(variance_recovery, file.path(paths$tables, "07_sim_lowrank_recovery_variances.csv"))
}

plot_status <- tibble(plot = character(), status = character())
plots_available <- setdiff(c("ggplot2"), rownames(installed.packages()))
if (length(plots_available) == 0L) {
  suppressPackageStartupMessages(library(ggplot2))

  empirical_long <- empirical_real |>
    pivot_longer(
      ends_with("_empirical_variance"),
      names_to = "component",
      values_to = "variance"
    ) |>
    mutate(component = sub("_empirical_variance", "", component))
  p_empirical <- ggplot(empirical_long, aes(reorder(observed_variable, variance), variance, fill = component)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "Empirical variance contribution", fill = "Component") +
    theme_minimal(base_size = 10)
  ggsave(file.path(paths$figures, "07_real_empirical_variance_decomposition.png"), p_empirical, width = 9, height = 10, dpi = 150)

  if (file.exists(file.path(paths$tables, "07_sim_diagonal_variance_vs_truth.csv"))) {
    sim_diag_recovery <- read_csv(file.path(paths$tables, "07_sim_diagonal_variance_vs_truth.csv"), show_col_types = FALSE)
    p_diag <- ggplot(sim_diag_recovery, aes(truth_value, mean, color = component)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.8) +
      labs(x = "True variance", y = "Posterior mean variance", color = "Component") +
      theme_minimal(base_size = 10)
    ggsave(file.path(paths$figures, "07_sim_diagonal_covariance_diag_vs_truth.png"), p_diag, width = 7, height = 5, dpi = 150)
  }

  if (file.exists(file.path(paths$tables, "07_sim_lowrank_recovery_loadings.csv"))) {
    loading_recovery <- read_csv(file.path(paths$tables, "07_sim_lowrank_recovery_loadings.csv"), show_col_types = FALSE)
    p_load <- ggplot(loading_recovery, aes(truth_value, mean, color = parameter)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.75) +
      labs(x = "True loading", y = "Posterior mean loading", color = "Matrix") +
      theme_minimal(base_size = 10)
    ggsave(file.path(paths$figures, "07_sim_lowrank_loadings_vs_truth.png"), p_load, width = 7, height = 5, dpi = 150)
  }

  if (file.exists(file.path(paths$tables, "07_sim_lowrank_recovery_variances.csv"))) {
    variance_recovery <- read_csv(file.path(paths$tables, "07_sim_lowrank_recovery_variances.csv"), show_col_types = FALSE)
    p_var <- ggplot(variance_recovery, aes(truth_value, mean, color = parameter)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.75) +
      labs(x = "True value", y = "Posterior mean", color = "Parameter") +
      theme_minimal(base_size = 10)
    ggsave(file.path(paths$figures, "07_sim_lowrank_variances_vs_truth.png"), p_var, width = 7, height = 5, dpi = 150)
  }

  plot_status <- tibble(
    plot = "diagnostic_plots",
    status = "summary_plots_written_trace_plots_skipped_to_avoid_large_draw_loading"
  )
} else {
  plot_status <- tibble(plot = "diagnostic_plots", status = "skipped_missing_ggplot2")
}

write_csv(plot_status, file.path(paths$diagnostics, "07_plot_status.csv"))
write_csv(truth_variance, file.path(paths$tables, "07_simulation_truth_variance_reference.csv"))

message("07_diagnostics_and_recovery_plots complete.")
