# Loading analysis, chain diagnostics, and recovery metrics for the
# low-rank-Sigma_a + diagonal-Sigma_b model (stages 09 and 10).

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/11_loading_analysis_and_chains.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "posterior", "tidyr"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
})

model_objects <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names  <- model_objects$variable_names
Q_a            <- simulation$rank_a

# ── Load posterior summaries ─────────────────────────────────────────────────

sim_summary_path  <- file.path(paths$tables, "09_sim_lowrank_a_diag_b_posterior_summary.csv")
real_summary_path <- file.path(paths$tables, "10_real_lowrank_a_diag_b_posterior_summary.csv")

sim_summary  <- read_csv_if_exists(sim_summary_path)
real_summary <- read_csv_if_exists(real_summary_path)

if (is.null(sim_summary) && is.null(real_summary)) {
  message("No posterior summaries found. Run stages 09 and/or 10 first.")
  quit(save = "no", status = 0)
}

chain_plot_dir <- file.path(paths$figures, "11_chain_plots")
dir.create(chain_plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper: save plot ─────────────────────────────────────────────────────────

save_plot <- function(p, filename, width = 12, height = 7) {
  ggplot2::ggsave(filename, p, width = width, height = height, dpi = 150)
  message("Saved: ", basename(filename))
}

# ── Posterior health summary helper ──────────────────────────────────────────

health_row <- function(summ, model_id) {
  if (is.null(summ)) return(tibble::tibble(model_id = model_id, max_rhat = NA_real_, min_ess_bulk = NA_real_, n_rhat_gt_105 = NA_integer_))
  tibble::tibble(
    model_id      = model_id,
    max_rhat      = max(summ$rhat, na.rm = TRUE),
    min_ess_bulk  = min(summ$ess_bulk, na.rm = TRUE),
    min_ess_tail  = min(summ$ess_tail, na.rm = TRUE),
    n_rhat_gt_101 = sum(!is.na(summ$rhat) & summ$rhat > 1.01),
    n_rhat_gt_105 = sum(!is.na(summ$rhat) & summ$rhat > 1.05)
  )
}

health <- dplyr::bind_rows(
  health_row(sim_summary,  "09_sim_lowrank_a_diag_b"),
  health_row(real_summary, "10_real_lowrank_a_diag_b")
)
readr::write_csv(health, file.path(paths$diagnostics, "11_posterior_health_new_model.csv"))
message("Posterior health:\n", paste(capture.output(print(health)), collapse = "\n"))

# ── Loading table function ────────────────────────────────────────────────────

make_and_save_loading_table <- function(summ, tag) {
  if (is.null(summ)) return(NULL)
  tbl <- build_loading_table(summ, feature_names, Q_a)
  out_path <- file.path(paths$tables, paste0(tag, "_lambda_a_loading_table.csv"))
  readr::write_csv(tbl, out_path)
  message("Loading table saved: ", basename(out_path))
  tbl
}

sim_loadings  <- make_and_save_loading_table(sim_summary,  "11_sim")
real_loadings <- make_and_save_loading_table(real_summary, "11_real")

# ── Heatmaps ─────────────────────────────────────────────────────────────────

if (!is.null(sim_loadings)) {
  p <- plot_loading_heatmap(sim_loadings, title = "Lambda_a loadings — simulated data (stage 09)")
  save_plot(p, file.path(paths$figures, "11_sim_lambda_a_heatmap.png"), width = 7, height = 14)
}

if (!is.null(real_loadings)) {
  p <- plot_loading_heatmap(real_loadings, title = "Lambda_a loadings — real La Liga data (stage 10)")
  save_plot(p, file.path(paths$figures, "11_real_lambda_a_heatmap.png"), width = 7, height = 14)
}

# ── Top-loading bar charts ────────────────────────────────────────────────────

if (!is.null(sim_loadings)) {
  p <- plot_top_loadings(sim_loadings, top_n = 12)
  save_plot(p, file.path(paths$figures, "11_sim_lambda_a_top_loadings.png"), width = 12, height = 8)
}

if (!is.null(real_loadings)) {
  p <- plot_top_loadings(real_loadings, top_n = 12)
  save_plot(p, file.path(paths$figures, "11_real_lambda_a_top_loadings.png"), width = 12, height = 8)
}

# ── Group radar / bar chart ───────────────────────────────────────────────────

if (!is.null(real_loadings)) {
  p <- plot_radar_by_group(real_loadings)
  save_plot(p, file.path(paths$figures, "11_real_lambda_a_group_radar.png"), width = 10, height = 6)
}

# ── Variance decomposition comparison ────────────────────────────────────────

make_variance_decomp <- function(summ, tag, feature_names) {
  if (is.null(summ)) return(NULL)
  sa <- extract_vector_summary(summ, "Sigma_a_diag", feature_names)
  vb <- extract_vector_summary(summ, "var_b",        feature_names)
  ve <- extract_vector_summary(summ, "var_e",        feature_names)
  pa <- extract_vector_summary(summ, "prop_a",       feature_names)
  pb <- extract_vector_summary(summ, "prop_b",       feature_names)
  pe <- extract_vector_summary(summ, "prop_e",       feature_names)

  if (is.null(sa) || is.null(pa)) return(NULL)

  tbl <- tibble::tibble(
    feature     = feature_names,
    sigma_a_diag_mean  = sa$mean,
    var_b_mean         = vb$mean,
    var_e_mean         = ve$mean,
    prop_a_mean        = pa$mean,
    prop_b_mean        = pb$mean,
    prop_e_mean        = pe$mean
  )
  out_path <- file.path(paths$tables, paste0(tag, "_variance_decomposition.csv"))
  readr::write_csv(tbl, out_path)
  message("Variance decomposition saved: ", basename(out_path))
  tbl
}

sim_vardecomp  <- make_variance_decomp(sim_summary,  "11_sim",  feature_names)
real_vardecomp <- make_variance_decomp(real_summary, "11_real", feature_names)

# ── Variance decomposition stacked bar (real data) ───────────────────────────

if (!is.null(real_vardecomp)) {
  groups <- feature_group_lookup()
  plot_df <- real_vardecomp |>
    dplyr::left_join(groups, by = "feature") |>
    dplyr::mutate(
      group = dplyr::coalesce(group, "Other"),
      group = factor(group, levels = group_order()),
      feature_label = gsub("per90_|rate_", "", feature),
      feature_label = gsub("_", " ", feature_label),
      feature_f = stats::reorder(feature_label, prop_a_mean)
    ) |>
    tidyr::pivot_longer(
      c(prop_a_mean, prop_b_mean, prop_e_mean),
      names_to = "component", values_to = "proportion"
    ) |>
    dplyr::mutate(component = dplyr::recode(component,
      prop_a_mean = "Player (Σ_a)", prop_b_mean = "Season (Σ_b)", prop_e_mean = "Residual"
    ))

  p_vd <- ggplot2::ggplot(plot_df,
    ggplot2::aes(x = proportion, y = feature_f, fill = component)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(values = c("Player (Σ_a)" = "#2166ac", "Season (Σ_b)" = "#f4a582", "Residual" = "#d6604d")) +
    ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::labs(
      title = "Variance decomposition — low-rank Σ_a model (real data)",
      x = "Proportion of total variance", y = NULL, fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      strip.text.y    = ggplot2::element_text(angle = 0, hjust = 0, size = 7, face = "bold"),
      legend.position = "top"
    )
  save_plot(p_vd, file.path(paths$figures, "11_real_variance_decomposition_lowrank.png"), width = 10, height = 16)
}

# ── Recovery metrics vs simulation truth ─────────────────────────────────────

truth_path <- file.path(paths$processed, "simulated_lowrank_truth.rds")

if (!is.null(sim_summary) && file.exists(truth_path)) {
  truth <- readRDS(truth_path)

  # Lambda_a recovery
  lambda_a_post <- extract_matrix_summary(sim_summary, "Lambda_a", feature_names)

  if (!is.null(lambda_a_post) && !is.null(truth$Lambda_a)) {
    true_lambda <- truth$Lambda_a
    P_t <- nrow(true_lambda)
    Q_t <- ncol(true_lambda)

    recovery_rows <- lapply(seq_len(Q_t), function(q) {
      post_q <- lambda_a_post |> dplyr::filter(factor_index == q)
      if (nrow(post_q) == 0L) return(NULL)
      true_q <- true_lambda[seq_len(nrow(post_q)), q]
      tibble::tibble(
        factor  = q,
        n       = length(true_q),
        bias    = mean(post_q$mean - true_q),
        mae     = mean(abs(post_q$mean - true_q)),
        rmse    = sqrt(mean((post_q$mean - true_q)^2)),
        corr    = cor(post_q$mean, true_q),
        cov90   = mean(post_q$q5 <= true_q & true_q <= post_q$q95)
      )
    })

    lambda_recovery <- dplyr::bind_rows(recovery_rows)
    message("Lambda_a recovery (stage 09 vs truth):\n",
            paste(capture.output(print(lambda_recovery)), collapse = "\n"))

    # Sigma_a_diag recovery
    sigma_a_post <- extract_vector_summary(sim_summary, "Sigma_a_diag", feature_names)
    true_sigma_a_diag <- diag(
      truth$Lambda_a %*% t(truth$Lambda_a) + diag(truth$psi_a_sd^2)
    )
    sigma_a_recovery <- tibble::tibble(
      parameter = "Sigma_a_diag",
      n    = length(true_sigma_a_diag),
      bias = mean(sigma_a_post$mean - true_sigma_a_diag),
      mae  = mean(abs(sigma_a_post$mean - true_sigma_a_diag)),
      rmse = sqrt(mean((sigma_a_post$mean - true_sigma_a_diag)^2)),
      corr = cor(sigma_a_post$mean, true_sigma_a_diag),
      cov90 = mean(sigma_a_post$q5 <= true_sigma_a_diag &
                   true_sigma_a_diag <= sigma_a_post$q95)
    )

    # sigma_b^2 vs truth Sigma_b diagonal
    vb_post <- extract_vector_summary(sim_summary, "var_b", feature_names)
    true_sigma_b_diag <- diag(
      truth$Lambda_b %*% t(truth$Lambda_b) + diag(truth$psi_b_sd^2)
    )
    sigma_b_recovery <- tibble::tibble(
      parameter = "var_b_vs_Sigma_b_diag",
      n    = length(true_sigma_b_diag),
      bias = mean(vb_post$mean - true_sigma_b_diag),
      mae  = mean(abs(vb_post$mean - true_sigma_b_diag)),
      rmse = sqrt(mean((vb_post$mean - true_sigma_b_diag)^2)),
      corr = cor(vb_post$mean, true_sigma_b_diag),
      cov90 = mean(vb_post$q5 <= true_sigma_b_diag &
                   true_sigma_b_diag <= vb_post$q95)
    )

    all_recovery <- dplyr::bind_rows(
      dplyr::mutate(lambda_recovery, parameter = paste0("Lambda_a_factor", factor)) |>
        dplyr::select(-factor),
      sigma_a_recovery,
      sigma_b_recovery
    )
    readr::write_csv(all_recovery, file.path(paths$tables, "11_sim_recovery_metrics.csv"))
    message("Recovery metrics saved.")

    # Scatter: posterior mean vs truth for Lambda_a
    scatter_df <- lambda_a_post |>
      dplyr::mutate(
        true_val = as.vector(true_lambda[observed_index, factor_index]),
        factor_label = paste0("Factor ", factor_index)
      )

    p_scatter <- ggplot2::ggplot(scatter_df,
      ggplot2::aes(x = true_val, y = mean, colour = factor_label)
    ) +
      ggplot2::geom_point(alpha = 0.7, size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = q5, ymax = q95), width = 0, alpha = 0.4
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
      ggplot2::facet_wrap(~factor_label, scales = "free") +
      ggplot2::labs(
        title  = "Lambda_a: posterior mean vs simulation truth (stage 09)",
        x      = "True loading", y = "Posterior mean ± 90% CI",
        colour = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")

    save_plot(p_scatter, file.path(paths$figures, "11_sim_lambda_a_vs_truth.png"), width = 10, height = 5)
  }
}

# ── Chain / trace plots ───────────────────────────────────────────────────────

make_chain_plots <- function(fit_path, model_id, tag, Q_a) {
  fit <- load_cmdstan_fit(fit_path = fit_path, model_id = model_id)
  if (is.null(fit)) {
    message("No fit available for ", model_id, " — skipping chain plots.")
    return(invisible(NULL))
  }

  sub_dir <- file.path(chain_plot_dir, tag)
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

  # Anchor loadings
  anchor_vars <- sprintf("lambda_a_diag[%d]", seq_len(Q_a))
  plot_trace_density_to_dir(fit, anchor_vars, file.path(sub_dir, "lambda_a_diag"))

  # Top free loadings by posterior mean absolute value (up to 8)
  summ <- fit$summary(variables = "lambda_a_free")
  if (!is.null(summ) && nrow(summ) > 0L) {
    top_free <- summ |>
      dplyr::arrange(dplyr::desc(abs(mean))) |>
      dplyr::slice_head(n = 8) |>
      dplyr::pull(variable)
    plot_trace_density_to_dir(fit, top_free, file.path(sub_dir, "lambda_a_free_top"))
  }

  # Uniqueness psi_a for first 6 features
  psi_vars <- sprintf("psi_a[%d]", seq_len(min(6L, length(feature_names))))
  plot_trace_density_to_dir(fit, psi_vars, file.path(sub_dir, "psi_a"))

  # sigma_b for 6 features with highest posterior mean season variance
  summ_sb <- fit$summary(variables = "sigma_b")
  if (!is.null(summ_sb) && nrow(summ_sb) >= 4L) {
    top_sb <- summ_sb |>
      dplyr::arrange(dplyr::desc(mean)) |>
      dplyr::slice_head(n = 6) |>
      dplyr::pull(variable)
    plot_trace_density_to_dir(fit, top_sb, file.path(sub_dir, "sigma_b_top"))
  }

  # sigma_e for first 6 features
  se_vars <- sprintf("sigma_e[%d]", seq_len(min(6L, length(feature_names))))
  plot_trace_density_to_dir(fit, se_vars, file.path(sub_dir, "sigma_e"))

  message("Chain plots saved to: ", sub_dir)
}

# Wrapper that saves to a specific directory prefix
plot_trace_density_to_dir <- function(fit, variables, prefix) {
  if (is.null(fit) || length(variables) == 0L) return(invisible(NULL))
  check_packages(c("posterior", "ggplot2", "tidyr"))
  tryCatch({
    draws <- posterior::as_draws_df(fit$draws(variables = variables))
    available <- intersect(variables, names(draws))
    if (length(available) == 0L) return(invisible(NULL))

    long <- draws |>
      dplyr::select(.chain, .iteration, dplyr::all_of(available)) |>
      tidyr::pivot_longer(dplyr::all_of(available), names_to = "parameter", values_to = "value")

    n_params <- length(available)
    ncols <- min(3L, n_params)
    nrows <- ceiling(n_params / ncols)

    trace_plot <- ggplot2::ggplot(long, ggplot2::aes(.iteration, value, color = factor(.chain))) +
      ggplot2::geom_line(linewidth = 0.25, alpha = 0.75) +
      ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = ncols) +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(x = "Iteration", y = "Value", color = "Chain") +
      ggplot2::theme_minimal(base_size = 9)

    density_plot <- ggplot2::ggplot(long, ggplot2::aes(value, fill = factor(.chain))) +
      ggplot2::geom_density(alpha = 0.35, colour = NA) +
      ggplot2::facet_wrap(~ parameter, scales = "free", ncol = ncols) +
      ggplot2::scale_fill_brewer(palette = "Set1") +
      ggplot2::labs(x = "Value", y = "Density", fill = "Chain") +
      ggplot2::theme_minimal(base_size = 9)

    height <- max(4, nrows * 2.5)
    ggplot2::ggsave(paste0(prefix, "_trace.png"),   trace_plot,   width = ncols * 3 + 1, height = height, dpi = 150)
    ggplot2::ggsave(paste0(prefix, "_density.png"), density_plot, width = ncols * 3 + 1, height = height, dpi = 150)
  }, error = function(e) {
    message("Chain plot failed for ", basename(prefix), ": ", conditionMessage(e))
  })
}

make_chain_plots(
  fit_path = file.path(paths$fits, "09_sim_lowrank_a_diag_b_fit.rds"),
  model_id = "09_sim_lowrank_a_diag_b",
  tag      = "09_sim",
  Q_a      = Q_a
)

make_chain_plots(
  fit_path = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
  model_id = "10_real_lowrank_a_diag_b",
  tag      = "10_real",
  Q_a      = Q_a
)

message("11_loading_analysis_and_chains complete.")
