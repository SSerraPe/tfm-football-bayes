# PCA-oriented post-processing of Lambda_a posterior draws.
#
# The lower-triangular Λ is an identification device, not an interpretable
# object.  Following the professor's postprocess_lowrank_diag.R, we
# eigendecompose B = ΛΛ' for each posterior draw and summarise the
# PCA-oriented loadings and player scores across draws.
#
# Reads:   fits/10_real_lowrank_a_diag_b_fit.rds  (and 09 sim if available)
# Writes:  outputs/tables/13_*.csv
#          outputs/figures/13_*.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/13_pca_postprocessing.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "tidyr"))

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
P              <- length(feature_names)
I_players      <- length(model_objects$player_levels)

N_USE <- 1000L   # number of posterior draws to use

# ── PCA helpers (from professor's postprocess_lowrank_diag.R) ─────────────────

pca_from_B <- function(B, K) {
  eg   <- eigen(B, symmetric = TRUE)
  vals <- pmax(eg$values[seq_len(K)], 0)
  U    <- eg$vectors[, seq_len(K), drop = FALSE]
  list(U = U, values = vals, loadings = U %*% diag(sqrt(vals), K))
}

fix_signs <- function(L) {
  for (h in seq_len(ncol(L))) {
    jmax <- which.max(abs(L[, h]))
    if (L[jmax, h] < 0) L[, h] <- -L[, h]
  }
  L
}

align_to_reference <- function(L, G, L_ref) {
  for (h in seq_len(ncol(L))) {
    if (sum(L[, h] * L_ref[, h]) < 0) {
      L[, h] <- -L[, h]
      G[, h] <- -G[, h]
    }
  }
  list(loadings = L, scores = G)
}

# ── Core post-processing loop ─────────────────────────────────────────────────

run_pca_postproc <- function(fit, tag, truth_Lambda = NULL) {
  if (is.null(fit)) {
    message("No fit for ", tag, " — skipping.")
    return(invisible(NULL))
  }

  # Variable names
  Lambda_vars <- as.vector(
    outer(seq_len(P), seq_len(Q_a),
          function(p, q) sprintf("Lambda_a[%d,%d]", p, q))
  )
  eta_vars <- as.vector(
    outer(seq_len(I_players), seq_len(Q_a),
          function(i, q) sprintf("eta_a_raw[%d,%d]", i, q))
  )

  message("Extracting Lambda_a draws for ", tag, " ...")
  lambda_mat <- tryCatch(
    posterior::as_draws_matrix(fit$draws(variables = "Lambda_a")),
    error = function(e) { message("  Lambda_a draw extraction failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(lambda_mat)) return(invisible(NULL))

  message("Extracting eta_a_raw draws for ", tag, " ...")
  eta_mat <- tryCatch(
    posterior::as_draws_matrix(fit$draws(variables = "eta_a_raw")),
    error = function(e) { message("  eta_a_raw draw extraction failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(eta_mat)) return(invisible(NULL))

  n_draws <- nrow(lambda_mat)
  set.seed(42L)
  use_idx <- sample(seq_len(n_draws), min(N_USE, n_draws))
  n_use   <- length(use_idx)

  # Step 1: build reference B from posterior mean
  message("Computing reference orientation from posterior mean B ...")
  B_ref <- matrix(0, P, P)
  for (s in use_idx) {
    L_s  <- matrix(lambda_mat[s, Lambda_vars], P, Q_a)
    B_ref <- B_ref + (L_s %*% t(L_s)) / n_use
  }
  pca_ref <- pca_from_B(B_ref, Q_a)
  L_ref   <- fix_signs(pca_ref$loadings)

  # Step 2: draw-wise PCA loadings and player scores
  message("Running draw-wise PCA post-processing (n = ", n_use, ") ...")
  L_pca_draws <- array(NA_real_, dim = c(n_use, P, Q_a))
  G_pca_draws <- array(NA_real_, dim = c(n_use, I_players, Q_a))
  eig_draws   <- matrix(NA_real_, n_use, Q_a)

  for (idx in seq_len(n_use)) {
    s    <- use_idx[idx]
    L_s  <- matrix(lambda_mat[s, Lambda_vars], P, Q_a)
    F_s  <- matrix(eta_mat[s, eta_vars],       I_players, Q_a)

    B_s   <- L_s %*% t(L_s)
    pca_s <- pca_from_B(B_s, Q_a)
    U_s   <- pca_s$U
    d_s   <- pca_s$values

    # Common player effect and PCA-oriented scores
    C_s <- F_s %*% t(L_s)
    safe_d <- pmax(d_s, 1e-10)
    G_s <- C_s %*% U_s %*% diag(1 / sqrt(safe_d), Q_a)

    aligned <- align_to_reference(pca_s$loadings, G_s, L_ref)
    L_pca_draws[idx, , ] <- aligned$loadings
    G_pca_draws[idx, , ] <- aligned$scores
    eig_draws[idx, ]     <- d_s
  }

  # Step 3: summarise loadings
  L_pca_mean  <- apply(L_pca_draws, c(2, 3), mean)
  L_pca_q05   <- apply(L_pca_draws, c(2, 3), quantile, probs = 0.05)
  L_pca_q95   <- apply(L_pca_draws, c(2, 3), quantile, probs = 0.95)

  # Step 4: summarise player scores
  G_pca_mean  <- apply(G_pca_draws, c(2, 3), mean)
  G_pca_sd    <- apply(G_pca_draws, c(2, 3), sd)
  G_pca_q05   <- apply(G_pca_draws, c(2, 3), quantile, probs = 0.05)
  G_pca_q95   <- apply(G_pca_draws, c(2, 3), quantile, probs = 0.95)

  # Step 5: proportion of common variance explained per PCA factor
  prop_common <- sweep(eig_draws, 1, rowSums(eig_draws), "/")

  # ── Save tables ─────────────────────────────────────────────────────────────

  groups <- feature_group_lookup()

  loading_rows <- lapply(seq_len(Q_a), function(q) {
    tibble::tibble(
      feature  = feature_names,
      factor   = paste0("PCA_Factor_", q),
      loading  = L_pca_mean[, q],
      q5       = L_pca_q05[, q],
      q95      = L_pca_q95[, q]
    )
  })
  loading_tbl <- dplyr::bind_rows(loading_rows) |>
    dplyr::left_join(groups, by = "feature") |>
    dplyr::mutate(
      group       = dplyr::coalesce(group, "Other"),
      group       = factor(group, levels = group_order()),
      abs_loading = abs(loading)
    )
  readr::write_csv(loading_tbl,
    file.path(paths$tables, paste0(tag, "_pca_loading_summary.csv")))
  message("PCA loading summary saved.")

  prop_tbl <- tibble::as_tibble(
    data.frame(
      factor   = paste0("PCA_Factor_", seq_len(Q_a)),
      median   = apply(prop_common, 2, median),
      q05      = apply(prop_common, 2, quantile, probs = 0.05),
      q95      = apply(prop_common, 2, quantile, probs = 0.95)
    )
  )
  readr::write_csv(prop_tbl,
    file.path(paths$tables, paste0(tag, "_pca_common_variance_prop.csv")))
  message("PCA common variance proportions saved.")

  score_tbl <- tibble::tibble(
    player_id = model_objects$player_levels
  )
  for (q in seq_len(Q_a)) {
    score_tbl[[paste0("pca_score_", q, "_mean")]] <- G_pca_mean[, q]
    score_tbl[[paste0("pca_score_", q, "_sd")]]   <- G_pca_sd[, q]
    score_tbl[[paste0("pca_score_", q, "_q05")]]  <- G_pca_q05[, q]
    score_tbl[[paste0("pca_score_", q, "_q95")]]  <- G_pca_q95[, q]
  }
  readr::write_csv(score_tbl,
    file.path(paths$tables, paste0(tag, "_pca_player_scores.csv")))
  message("PCA player scores saved.")

  # ── Figures ──────────────────────────────────────────────────────────────────

  save_plot <- function(p, fname, w = 10, h = 7) {
    ggplot2::ggsave(fname, p, width = w, height = h, dpi = 150)
    message("Saved: ", basename(fname))
  }

  # Heatmap of PCA loadings
  heat_df <- loading_tbl |>
    dplyr::arrange(group, feature) |>
    dplyr::mutate(
      feature_label = gsub("per90_|rate_", "", feature),
      feature_label = gsub("_", " ", feature_label),
      feature_f     = factor(feature_label, levels = rev(unique(feature_label)))
    )
  p_heat <- ggplot2::ggplot(heat_df,
    ggplot2::aes(x = factor, y = feature_f, fill = loading)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", loading)),
      size = 2.3, colour = "white", fontface = "bold") +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d",
      midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish, name = "Loading") +
    ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(title = paste0("PCA-oriented loadings — ", tag),
      x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, size = 8, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 7),
      panel.grid   = ggplot2::element_blank()
    )
  save_plot(p_heat,
    file.path(paths$figures, paste0(tag, "_pca_loading_heatmap.png")), w = 7, h = 14)

  # CI dot plot for PCA loadings
  ci_df <- loading_tbl |>
    dplyr::mutate(
      feature_label = gsub("per90_|rate_", "", feature),
      feature_label = gsub("_", " ", feature_label),
      feature_f     = stats::reorder(feature_label, loading)
    )
  p_ci <- ggplot2::ggplot(ci_df,
    ggplot2::aes(x = loading, y = feature_f, colour = group)) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = q5, xmax = q95),
      height = 0.25, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = group_colours(), name = "Group") +
    ggplot2::facet_wrap(~factor) +
    ggplot2::labs(
      title = paste0("PCA-oriented loadings with 90% CI — ", tag),
      x = "PCA loading", y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "bottom",
      axis.text.y = ggplot2::element_text(size = 7))
  save_plot(p_ci,
    file.path(paths$figures, paste0(tag, "_pca_loading_ci_plot.png")), w = 14, h = 12)

  # Boxplot of proportion of common variance
  prop_long <- as.data.frame(prop_common)
  colnames(prop_long) <- paste0("PCA_Factor_", seq_len(Q_a))
  prop_long <- tidyr::pivot_longer(prop_long, dplyr::everything(),
    names_to = "factor", values_to = "proportion")
  p_prop <- ggplot2::ggplot(prop_long,
    ggplot2::aes(x = factor, y = proportion, fill = factor)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
    ggplot2::labs(
      title = paste0("Common player variance explained per PCA factor — ", tag),
      x = "PCA factor", y = "Proportion of common variance") +
    ggplot2::theme_minimal(base_size = 11)
  save_plot(p_prop,
    file.path(paths$figures, paste0(tag, "_pca_common_variance_boxplot.png")), w = 6, h = 5)

  # Player score scatter (Factor 1 vs Factor 2)
  if (Q_a >= 2L) {
    player_map_path <- file.path(paths$processed, "player_map.csv")
    score_plot_df   <- score_tbl
    if (file.exists(player_map_path)) {
      player_map    <- readr::read_csv(player_map_path, show_col_types = FALSE)
      score_plot_df <- dplyr::left_join(score_plot_df, player_map, by = "player_id")
    }

    top_f1  <- score_plot_df |> dplyr::slice_max(pca_score_1_mean, n = 12)
    top_f2  <- score_plot_df |> dplyr::slice_max(pca_score_2_mean, n = 12)
    bot_f1  <- score_plot_df |> dplyr::slice_min(pca_score_1_mean, n = 8)
    bot_f2  <- score_plot_df |> dplyr::slice_min(pca_score_2_mean, n = 8)
    to_label <- dplyr::bind_rows(top_f1, top_f2, bot_f1, bot_f2) |> dplyr::distinct()

    player_col <- if ("player_name" %in% names(score_plot_df)) "player_name" else "player_id"

    p_scores <- ggplot2::ggplot(score_plot_df,
      ggplot2::aes(x = pca_score_1_mean, y = pca_score_2_mean)) +
      ggplot2::geom_point(alpha = 0.25, size = 1, colour = "#2166ac") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_point(data = to_label, colour = "#d6604d", size = 2.5) +
      ggplot2::labs(
        title = paste0("PCA-oriented player scores — ", tag),
        x = "PCA Factor 1 score",
        y = "PCA Factor 2 score") +
      ggplot2::theme_minimal(base_size = 10)

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p_scores <- p_scores +
        ggrepel::geom_text_repel(
          data = to_label,
          ggplot2::aes(label = .data[[player_col]]),
          size = 2.5, max.overlaps = 20)
    }
    save_plot(p_scores,
      file.path(paths$figures, paste0(tag, "_pca_player_scores_scatter.png")), w = 10, h = 8)
  }

  # ── Recovery validation for simulated data ───────────────────────────────────

  if (!is.null(truth_Lambda)) {
    B_true     <- truth_Lambda %*% t(truth_Lambda)
    pca_true   <- pca_from_B(B_true, Q_a)
    L_pca_true <- fix_signs(pca_true$loadings)
    # Align truth to reference
    dummy_G    <- matrix(0, I_players, Q_a)
    aligned_t  <- align_to_reference(L_pca_true, dummy_G, L_ref)
    L_pca_true <- aligned_t$loadings

    recovery <- tibble::tibble(
      factor         = paste0("PCA_Factor_", seq_len(Q_a)),
      cor_loading    = vapply(seq_len(Q_a), function(q)
        cor(L_pca_true[, q], L_pca_mean[, q]), numeric(1)),
      rmse_loading   = vapply(seq_len(Q_a), function(q)
        sqrt(mean((L_pca_mean[, q] - L_pca_true[, q])^2)), numeric(1)),
      prop_cv_truth  = pca_true$values / sum(pca_true$values),
      prop_cv_median = apply(prop_common, 2, median)
    )
    readr::write_csv(recovery,
      file.path(paths$tables, paste0(tag, "_pca_recovery.csv")))
    message("PCA recovery metrics saved:\n",
      paste(capture.output(print(recovery)), collapse = "\n"))

    # Scatter: true vs posterior mean PCA loadings
    scatter_df <- dplyr::bind_rows(lapply(seq_len(Q_a), function(q) {
      tibble::tibble(
        factor    = paste0("PCA_Factor_", q),
        true_val  = L_pca_true[, q],
        post_mean = L_pca_mean[, q],
        post_q05  = L_pca_q05[, q],
        post_q95  = L_pca_q95[, q]
      )
    }))
    p_rec <- ggplot2::ggplot(scatter_df,
      ggplot2::aes(x = true_val, y = post_mean, colour = factor)) +
      ggplot2::geom_point(alpha = 0.7, size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = post_q05, ymax = post_q95), width = 0, alpha = 0.4) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::facet_wrap(~factor, scales = "free") +
      ggplot2::labs(
        title = paste0("PCA loading recovery — ", tag),
        x = "True PCA loading", y = "Posterior mean ± 90% CI",
        colour = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "none")
    save_plot(p_rec,
      file.path(paths$figures, paste0(tag, "_pca_loading_recovery.png")), w = 10, h = 5)
  }

  message(tag, " PCA post-processing complete.")
  invisible(list(L_pca_mean = L_pca_mean, G_pca_mean = G_pca_mean,
                 eig_draws = eig_draws))
}

# ── Run on real data (stage 10) ───────────────────────────────────────────────

message("=== PCA post-processing: real data (stage 10) ===")
real_fit <- load_cmdstan_fit(
  fit_path = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
  model_id = "10_real_lowrank_a_diag_b"
)
run_pca_postproc(real_fit, tag = "13_real")

# ── Run on simulated data (stage 09), with recovery validation ────────────────

message("=== PCA post-processing: simulated data (stage 09) ===")
sim_fit <- load_cmdstan_fit(
  fit_path = file.path(paths$fits, "09_sim_lowrank_a_diag_b_fit.rds"),
  model_id = "09_sim_lowrank_a_diag_b"
)

truth_path <- file.path(paths$processed, "simulated_lowrank_truth.rds")
truth_Lambda <- if (file.exists(truth_path)) {
  truth <- readRDS(truth_path)
  if (!is.null(truth$Lambda_a)) truth$Lambda_a else NULL
} else NULL

run_pca_postproc(sim_fit, tag = "13_sim", truth_Lambda = truth_Lambda)

message("13_pca_postprocessing complete.")
