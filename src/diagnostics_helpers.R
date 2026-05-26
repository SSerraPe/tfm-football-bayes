# Helpers for posterior diagnostics and recovery summaries.

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

read_fit_if_exists <- function(path) {
  model_id <- sub("_fit[.]rds$", "", basename(path))
  load_cmdstan_fit(fit_path = path, model_id = model_id)
}

extract_vector_summary <- function(summary_table, parameter, labels) {
  if (is.null(summary_table)) return(NULL)
  prefix <- paste0(parameter, "[")
  summary_table |>
    dplyr::filter(startsWith(variable, prefix)) |>
    dplyr::mutate(
      index = as.integer(sub("\\].*", "", sub(prefix, "", variable, fixed = TRUE))),
      observed_variable = labels[index],
      parameter = parameter
    ) |>
    dplyr::select(parameter, index, observed_variable, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)
}

extract_matrix_summary <- function(summary_table, parameter, labels) {
  if (is.null(summary_table)) return(NULL)
  prefix <- paste0(parameter, "[")
  summary_table |>
    dplyr::filter(startsWith(variable, prefix)) |>
    dplyr::mutate(
      index_text = substring(variable, nchar(prefix) + 1, nchar(variable) - 1),
      observed_index = as.integer(sub(",.*", "", index_text)),
      factor_index = as.integer(sub(".*,", "", index_text)),
      observed_variable = labels[observed_index],
      parameter = parameter
    ) |>
    dplyr::select(
      parameter, observed_index, observed_variable, factor_index,
      mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail
    )
}

empirical_variance_decomposition <- function(model_objects) {
  Y <- model_objects$Y
  player_index <- model_objects$player_index
  season_index <- model_objects$season_index
  variable_names <- model_objects$variable_names

  rows <- vector("list", ncol(Y))
  for (p in seq_len(ncol(Y))) {
    y <- Y[, p]
    player_mean <- ave(y, player_index, FUN = mean)
    season_mean <- ave(y, season_index, FUN = mean)
    fitted <- player_mean + season_mean - mean(y)
    residual <- y - fitted
    rows[[p]] <- tibble::tibble(
      observed_variable = variable_names[p],
      player_empirical_variance = stats::var(player_mean),
      season_empirical_variance = stats::var(season_mean),
      residual_empirical_variance = stats::var(residual)
    )
  }
  dplyr::bind_rows(rows)
}

plot_trace_density <- function(fit, variables, prefix) {
  if (is.null(fit) || length(variables) == 0L) return(invisible(NULL))
  check_packages(c("posterior", "ggplot2", "tidyr"))
  draws <- posterior::as_draws_df(fit$draws(variables = variables))
  available <- intersect(variables, names(draws))
  if (length(available) == 0L) return(invisible(NULL))

  long <- draws |>
    dplyr::select(.chain, .iteration, dplyr::all_of(available)) |>
    tidyr::pivot_longer(dplyr::all_of(available), names_to = "parameter", values_to = "value")

  trace_plot <- ggplot2::ggplot(long, ggplot2::aes(.iteration, value, color = factor(.chain))) +
    ggplot2::geom_line(linewidth = 0.25, alpha = 0.75) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = "Value", color = "Chain") +
    ggplot2::theme_minimal(base_size = 10)

  density_plot <- ggplot2::ggplot(long, ggplot2::aes(value, fill = factor(.chain))) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::facet_wrap(~ parameter, scales = "free") +
    ggplot2::labs(x = "Value", y = "Density", fill = "Chain") +
    ggplot2::theme_minimal(base_size = 10)

  ggplot2::ggsave(file.path(paths$figures, paste0(prefix, "_trace.png")), trace_plot, width = 10, height = 6, dpi = 150)
  ggplot2::ggsave(file.path(paths$figures, paste0(prefix, "_density.png")), density_plot, width = 10, height = 6, dpi = 150)
  invisible(NULL)
}
