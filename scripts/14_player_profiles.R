# Player profile analysis: posterior means of a_i, clustered heatmap, PCA scatter.
#
# Uses the stage 10 real-data fit.  The full player effect is
#   a_i = Λ (η_i - η̄) + ψ ⊙ (z_{u,i} - z̄_u)
# For visualisation we use the low-rank component (common part = Λ η_i)
# which captures the shared factor structure and is compact (I × Q_a).
# For the clustered heatmap we reconstruct the full a_i in R^P.
#
# PCA-oriented scores from stage 13 are used for the scatter plot; they are
# recomputed inline if the stage-13 CSV is not available.
#
# Reads:   fits/10_real_lowrank_a_diag_b_fit.rds
#          outputs/tables/13_real_pca_player_scores.csv  (if available)
# Writes:  outputs/tables/14_player_effect_means.csv
#          outputs/tables/14_player_pca_scores.csv       (if stage 13 missing)
#          outputs/figures/14_player_profile_heatmap.png
#          outputs/figures/14_player_pca_scatter.png
#          outputs/figures/14_player_score_uncertainty.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/14_player_profiles.R"
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

N_USE <- 500L   # posterior draws for player effect reconstruction (expensive)

save_plot <- function(p, fname, w = 10, h = 7) {
  ggplot2::ggsave(fname, p, width = w, height = h, dpi = 150)
  message("Saved: ", basename(fname))
}

# ── Load fit ──────────────────────────────────────────────────────────────────

fit <- load_cmdstan_fit(
  fit_path = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
  model_id = "10_real_lowrank_a_diag_b"
)
if (is.null(fit)) {
  message("No fit available for stage 10 — skipping 14_player_profiles.")
  quit(save = "no", status = 0)
}

# ── Extract draws ─────────────────────────────────────────────────────────────

message("Extracting parameter draws ...")
Lambda_vars <- as.vector(
  outer(seq_len(P), seq_len(Q_a),
        function(p, q) sprintf("Lambda_a[%d,%d]", p, q))
)
eta_vars <- as.vector(
  outer(seq_len(I_players), seq_len(Q_a),
        function(i, q) sprintf("eta_a_raw[%d,%d]", i, q))
)
z_a_vars <- as.vector(
  outer(seq_len(I_players), seq_len(P),
        function(i, p) sprintf("z_a_raw[%d,%d]", i, p))
)
psi_a_vars  <- sprintf("psi_a[%d]", seq_len(P))

lambda_mat <- posterior::as_draws_matrix(fit$draws(variables = "Lambda_a"))
eta_mat    <- posterior::as_draws_matrix(fit$draws(variables = "eta_a_raw"))
z_a_mat    <- posterior::as_draws_matrix(fit$draws(variables = "z_a_raw"))
psi_a_mat  <- posterior::as_draws_matrix(fit$draws(variables = "psi_a"))

n_draws <- nrow(lambda_mat)
set.seed(42L)
use_idx <- sample(seq_len(n_draws), min(N_USE, n_draws))
n_use   <- length(use_idx)
message("Using ", n_use, " posterior draws for player effect reconstruction.")

# ── Reconstruct full player effects a_i ───────────────────────────────────────
# a_i[p] = sum_q Lambda_a[p,q] * (eta_a_raw[i,q] - mean_i(eta_a_raw[i,q]))
#         + psi_a[p] * (z_a_raw[i,p] - mean_i(z_a_raw[i,p]))

message("Reconstructing posterior mean player effects (", I_players, " players × ", P, " features) ...")
A_mean <- matrix(0, I_players, P)
A_sq   <- matrix(0, I_players, P)

for (idx in seq_len(n_use)) {
  s       <- use_idx[idx]
  L_s     <- matrix(lambda_mat[s, Lambda_vars], P, Q_a)
  eta_s   <- matrix(eta_mat[s, eta_vars],       I_players, Q_a)
  z_a_s   <- matrix(z_a_mat[s, z_a_vars],       I_players, P)
  psi_a_s <- as.numeric(psi_a_mat[s, psi_a_vars])

  # Column-centre (mirrors Stan's transformed parameters)
  eta_c <- sweep(eta_s, 2, colMeans(eta_s), "-")
  z_a_c <- sweep(z_a_s, 2, colMeans(z_a_s), "-")

  A_s <- eta_c %*% t(L_s) + sweep(z_a_c, 2, psi_a_s, "*")
  A_mean <- A_mean + A_s / n_use
  A_sq   <- A_sq   + A_s^2 / n_use
}

A_sd <- sqrt(pmax(A_sq - A_mean^2, 0))

# Save full player effect table
player_tbl <- as.data.frame(A_mean)
colnames(player_tbl) <- feature_names
player_tbl <- tibble::add_column(
  tibble::as_tibble(player_tbl),
  player_id = model_objects$player_levels,
  .before = 1
)

# Optionally join player names
player_map_path <- file.path(paths$processed, "player_map.csv")
if (file.exists(player_map_path)) {
  player_map <- readr::read_csv(player_map_path, show_col_types = FALSE)
  player_tbl <- dplyr::left_join(player_tbl, player_map, by = "player_id")
}
readr::write_csv(player_tbl,
  file.path(paths$tables, "14_player_effect_means.csv"))
message("Player effect means saved.")

# ── PCA-oriented scores ───────────────────────────────────────────────────────

scores_path <- file.path(paths$tables, "13_real_pca_player_scores.csv")
if (file.exists(scores_path)) {
  pca_scores <- readr::read_csv(scores_path, show_col_types = FALSE)
  message("Loaded PCA scores from stage 13.")
} else {
  message("Stage 13 PCA scores not found — computing inline (reduced draws) ...")

  pca_from_B <- function(B, K) {
    eg   <- eigen(B, symmetric = TRUE)
    vals <- pmax(eg$values[seq_len(K)], 0)
    U    <- eg$vectors[, seq_len(K), drop = FALSE]
    list(U = U, values = vals)
  }
  fix_signs <- function(L) {
    for (h in seq_len(ncol(L))) {
      jmax <- which.max(abs(L[, h]))
      if (L[jmax, h] < 0) L[, h] <- -L[, h]
    }
    L
  }

  # Reference orientation
  B_ref <- matrix(0, P, P)
  for (idx in seq_len(n_use)) {
    s   <- use_idx[idx]
    L_s <- matrix(lambda_mat[s, Lambda_vars], P, Q_a)
    B_ref <- B_ref + (L_s %*% t(L_s)) / n_use
  }
  pca_ref <- pca_from_B(B_ref, Q_a)
  U_ref   <- fix_signs(pca_ref$U %*% diag(sqrt(pca_ref$values), Q_a))

  G_pca_sum  <- matrix(0, I_players, Q_a)
  G_pca_sum2 <- matrix(0, I_players, Q_a)

  for (idx in seq_len(n_use)) {
    s    <- use_idx[idx]
    L_s  <- matrix(lambda_mat[s, Lambda_vars], P, Q_a)
    eta_s <- matrix(eta_mat[s, eta_vars],      I_players, Q_a)

    eta_c <- sweep(eta_s, 2, colMeans(eta_s), "-")
    C_s   <- eta_c %*% t(L_s)

    B_s   <- L_s %*% t(L_s)
    pca_s <- pca_from_B(B_s, Q_a)
    U_s   <- pca_s$U
    d_s   <- pmax(pca_s$values, 1e-10)

    G_s <- C_s %*% U_s %*% diag(1 / sqrt(d_s), Q_a)
    # Align signs to reference
    for (h in seq_len(Q_a)) {
      if (sum(G_s[, h] * U_ref[, h]) < 0) G_s[, h] <- -G_s[, h]
    }

    G_pca_sum  <- G_pca_sum  + G_s / n_use
    G_pca_sum2 <- G_pca_sum2 + G_s^2 / n_use
  }

  G_mean <- G_pca_sum
  G_sd   <- sqrt(pmax(G_pca_sum2 - G_mean^2, 0))

  pca_scores <- tibble::tibble(player_id = model_objects$player_levels)
  for (q in seq_len(Q_a)) {
    pca_scores[[paste0("pca_score_", q, "_mean")]] <- G_mean[, q]
    pca_scores[[paste0("pca_score_", q, "_sd")]]   <- G_sd[, q]
  }
  readr::write_csv(pca_scores,
    file.path(paths$tables, "14_player_pca_scores.csv"))
  message("Inline PCA scores saved.")
}

# Merge player names if available
if (file.exists(player_map_path)) {
  player_map  <- readr::read_csv(player_map_path, show_col_types = FALSE)
  pca_scores  <- dplyr::left_join(pca_scores, player_map, by = "player_id")
}
player_col <- if ("player_name" %in% names(pca_scores)) "player_name" else "player_id"

# ── Clustered heatmap (pheatmap or ggplot fallback) ───────────────────────────

# Select players to display: top 200 by L2 norm of player effect
norms <- sqrt(rowSums(A_mean^2))
top_idx <- order(norms, decreasing = TRUE)[seq_len(min(200L, I_players))]

A_sub <- A_mean[top_idx, ]
rownames(A_sub) <- model_objects$player_levels[top_idx]
colnames(A_sub) <- feature_names

# Feature annotation (group)
groups <- feature_group_lookup()
feat_ann <- data.frame(Group = groups$group[match(feature_names, groups$feature)],
                       row.names = feature_names)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  group_pal <- group_colours()
  ann_colors <- list(Group = group_pal[intersect(names(group_pal), unique(feat_ann$Group))])

  png(file.path(paths$figures, "14_player_profile_heatmap.png"),
      width = 1600, height = 2000, res = 130)
  pheatmap::pheatmap(
    t(A_sub),
    cluster_rows   = TRUE,
    cluster_cols   = TRUE,
    show_rownames  = TRUE,
    show_colnames  = FALSE,
    annotation_row = feat_ann,
    annotation_colors = ann_colors,
    color          = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
    fontsize_row   = 7,
    main           = "Posterior mean player effects (a_i) — top 200 players by norm",
    border_color   = NA
  )
  dev.off()
  message("Saved: 14_player_profile_heatmap.png")
} else {
  message("pheatmap not available — producing ggplot heatmap (first 50 players).")
  heat_df <- as.data.frame(A_sub[seq_len(min(50L, nrow(A_sub))), ]) |>
    tibble::rownames_to_column("player") |>
    tidyr::pivot_longer(-player, names_to = "feature", values_to = "effect") |>
    dplyr::left_join(groups, by = "feature") |>
    dplyr::mutate(group = dplyr::coalesce(group, "Other"),
                  feature_f = factor(feature,
                    levels = groups$feature[groups$feature %in% feature_names]))

  p_heat <- ggplot2::ggplot(heat_df,
    ggplot2::aes(x = player, y = feature_f, fill = effect)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0) +
    ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(title = "Posterior mean player effects (top 50)", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 4),
      strip.text.y = ggplot2::element_text(angle = 0, size = 6))
  save_plot(p_heat,
    file.path(paths$figures, "14_player_profile_heatmap.png"), w = 20, h = 14)
}

# ── Player PCA scatter ────────────────────────────────────────────────────────

if (Q_a >= 2L &&
    all(c("pca_score_1_mean", "pca_score_2_mean") %in% names(pca_scores))) {
  top_f1   <- pca_scores |> dplyr::slice_max(pca_score_1_mean, n = 12)
  top_f2   <- pca_scores |> dplyr::slice_max(pca_score_2_mean, n = 12)
  bot_f1   <- pca_scores |> dplyr::slice_min(pca_score_1_mean, n = 8)
  bot_f2   <- pca_scores |> dplyr::slice_min(pca_score_2_mean, n = 8)
  to_label <- dplyr::bind_rows(top_f1, top_f2, bot_f1, bot_f2) |> dplyr::distinct()

  p_scatter <- ggplot2::ggplot(pca_scores,
    ggplot2::aes(x = pca_score_1_mean, y = pca_score_2_mean)) +
    ggplot2::geom_point(alpha = 0.25, size = 1, colour = "#2166ac") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(data = to_label, colour = "#d6604d", size = 2.5) +
    ggplot2::labs(
      title = "Player archetypes in PCA factor space (posterior mean scores)",
      x = "PCA Factor 1 score", y = "PCA Factor 2 score") +
    ggplot2::theme_minimal(base_size = 10)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_scatter <- p_scatter +
      ggrepel::geom_text_repel(
        data = to_label,
        ggplot2::aes(label = .data[[player_col]]),
        size = 2.5, max.overlaps = 20)
  }
  save_plot(p_scatter,
    file.path(paths$figures, "14_player_pca_scatter.png"), w = 10, h = 8)
}

# ── Player score uncertainty ──────────────────────────────────────────────────

if (Q_a >= 1L && "pca_score_1_sd" %in% names(pca_scores)) {
  # Scatter: mean score vs posterior SD — shows which players are uncertain
  unc_df <- pca_scores |>
    dplyr::select(dplyr::starts_with("pca_score_")) |>
    tidyr::pivot_longer(
      dplyr::everything(),
      names_to  = c("factor", ".value"),
      names_pattern = "pca_score_(\\d+)_(.*)"
    ) |>
    dplyr::mutate(factor = paste0("PCA Factor ", factor))

  if (all(c("mean", "sd") %in% names(unc_df))) {
    p_unc <- ggplot2::ggplot(unc_df,
      ggplot2::aes(x = abs(mean), y = sd)) +
      ggplot2::geom_point(alpha = 0.2, size = 0.8, colour = "#2166ac") +
      ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "#d6604d",
        linewidth = 0.8) +
      ggplot2::facet_wrap(~factor) +
      ggplot2::labs(
        title = "Player score uncertainty vs magnitude",
        x = "|Posterior mean score|",
        y = "Posterior SD") +
      ggplot2::theme_minimal(base_size = 10)
    save_plot(p_unc,
      file.path(paths$figures, "14_player_score_uncertainty.png"), w = 10, h = 5)
  }
}

message("14_player_profiles complete.")
