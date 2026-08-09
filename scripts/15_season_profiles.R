# Season profile analysis: temporal evolution of b_j with credible intervals.
#
# Computes posterior mean and 90% CI of season effects b_j for each feature,
# then plots temporal line charts for features with the highest ICC_season.
#
# b_j[p] = sigma_b[p] * (z_b[j, p] - z_b_bar[p])   (mirrors Stan centering)
#
# Reads:   fits/10_real_lowrank_a_diag_b_fit.rds
#          outputs/tables/11_real_variance_decomposition.csv  (ICC_season ranking)
# Writes:  outputs/tables/15_season_effect_means.csv
#          outputs/figures/15_season_profiles_grid.png
#          outputs/figures/15_season_profile_heatmap.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/15_season_profiles.R"
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
season_levels  <- model_objects$season_levels
P              <- length(feature_names)
S              <- length(season_levels)

# Back-transformation helper ---------------------------------------------------
# Converts a season effect b_j (Z-score space) to a deviation in original units.
# center/scale are the post-transform, pre-Z-scale mean and SD (from scaling_parameters).
back_transform_season_effect <- function(b, center, scale, type) {
  y_t <- b * scale + center
  avg <- ifelse(type == "log1p", expm1(center),
         ifelse(type == "sqrt",  pmax(center, 0)^2, center))
  x   <- ifelse(type == "log1p", expm1(y_t),
         ifelse(type == "sqrt",  pmax(y_t, 0)^2,    y_t))
  x - avg
}

sp <- if (!is.null(model_objects$scaling_parameters) &&
          "transform_type" %in% names(model_objects$scaling_parameters)) {
  model_objects$scaling_parameters
} else {
  sp_path <- file.path(paths$processed, "scaling_parameters.csv")
  if (file.exists(sp_path)) {
    readr::read_csv(sp_path, show_col_types = FALSE)
  } else {
    NULL
  }
}

N_USE <- 1000L

save_plot <- function(p, fname, w = 12, h = 8) {
  ggplot2::ggsave(fname, p, width = w, height = h, dpi = 150)
  message("Saved: ", basename(fname))
}

# ── Load fit ──────────────────────────────────────────────────────────────────

# Fit selection: prefer the t-model at K* (from stage 32 ELPD, or BFA_SIM_RANK_A),
# then the stage-28 production fit, then the stage-10 Normal baseline.
rank_a_env <- Sys.getenv("BFA_SIM_RANK_A", unset = "")
chosen_k <- if (nchar(rank_a_env) > 0L) {
  as.integer(rank_a_env)
} else {
  elpd_path <- file.path(paths$tables, "32_player_holdout_elpd.csv")
  if (file.exists(elpd_path)) {
    et <- readr::read_csv(elpd_path, show_col_types = FALSE) |>
      dplyr::filter(!is.na(elpd_mv))
    if (nrow(et) > 0L) et$K[which.max(et$elpd_mv)] else NA_integer_
  } else NA_integer_
}

candidate_ids <- if (!is.na(chosen_k)) {
  c(sprintf("28_real_lowrank_a_diag_b_t_k%d", chosen_k),
    sprintf("34a_real_lowrank_a_diag_b_t_mv_k%d", chosen_k),
    "28_real_lowrank_a_diag_b_t_k3",
    "10_real_lowrank_a_diag_b")
} else {
  c("28_real_lowrank_a_diag_b_t_k3", "10_real_lowrank_a_diag_b")
}

fit <- NULL
for (.id in candidate_ids) {
  fit <- load_cmdstan_fit(
    fit_path = file.path(paths$fits, paste0(.id, "_fit.rds")),
    model_id = .id
  )
  if (!is.null(fit)) { message("Loaded fit: ", .id); break }
}
rm(.id)
if (is.null(fit)) {
  message("No fit available — skipping 15_season_profiles. Run stage 28 or set BFA_SIM_RANK_A.")
  quit(save = "no", status = 0)
}

# ── Extract draws ─────────────────────────────────────────────────────────────

message("Extracting sigma_b and z_b draws ...")
sigma_b_mat <- posterior::as_draws_matrix(fit$draws(variables = "sigma_b"))
z_b_mat     <- posterior::as_draws_matrix(fit$draws(variables = "z_b"))

sigma_b_vars <- sprintf("sigma_b[%d]", seq_len(P))
z_b_vars     <- as.vector(
  outer(seq_len(S), seq_len(P),
        function(j, p) sprintf("z_b[%d,%d]", j, p))
)

n_draws <- nrow(sigma_b_mat)
set.seed(42L)
use_idx <- sample(seq_len(n_draws), min(N_USE, n_draws))
n_use   <- length(use_idx)
message("Using ", n_use, " posterior draws for season effect reconstruction.")

# ── Reconstruct b_j: posterior mean and 90% CI ────────────────────────────────
# Accumulate E[b_j | data] across draws

message("Reconstructing season effects b_j (", S, " seasons × ", P, " features) ...")
B_sum  <- matrix(0, S, P)
B_sum2 <- matrix(0, S, P)
B_q05  <- array(NA_real_, dim = c(n_use, S, P))

for (idx in seq_len(n_use)) {
  s         <- use_idx[idx]
  sigma_b_s <- as.numeric(sigma_b_mat[s, sigma_b_vars])
  z_b_s     <- matrix(z_b_mat[s, z_b_vars], S, P)

  # Column-centre to mirror Stan's z_b_bar centering
  z_b_c <- sweep(z_b_s, 2, colMeans(z_b_s), "-")
  B_s   <- sweep(z_b_c, 2, sigma_b_s, "*")

  B_sum  <- B_sum  + B_s / n_use
  B_sum2 <- B_sum2 + B_s^2 / n_use
  B_q05[idx, , ] <- B_s
}

B_mean <- B_sum
B_sd   <- sqrt(pmax(B_sum2 - B_mean^2, 0))

# Compute quantiles from stored draws
B_lo <- apply(B_q05, c(2, 3), quantile, probs = 0.05)
B_hi <- apply(B_q05, c(2, 3), quantile, probs = 0.95)

# ── Save table (Z-score scale + original-unit deviations) ─────────────────────

season_rows <- lapply(seq_len(S), function(j) {
  tibble::tibble(
    season  = season_levels[j],
    feature = feature_names,
    b_mean  = B_mean[j, ],
    b_sd    = B_sd[j, ],
    b_q05   = B_lo[j, ],
    b_q95   = B_hi[j, ]
  )
})
season_tbl <- dplyr::bind_rows(season_rows)

# Add original-unit deviations if transform info is available
if (!is.null(sp) && "transform_type" %in% names(sp)) {
  sp_lookup <- dplyr::rename(sp, feature = variable)
  season_tbl <- season_tbl |>
    dplyr::left_join(sp_lookup, by = "feature") |>
    dplyr::mutate(
      b_orig_mean = back_transform_season_effect(b_mean, center, scale, transform_type),
      b_orig_q05  = back_transform_season_effect(b_q05,  center, scale, transform_type),
      b_orig_q95  = back_transform_season_effect(b_q95,  center, scale, transform_type)
    ) |>
    dplyr::select(-center, -scale, -transform_type)
}

readr::write_csv(season_tbl,
  file.path(paths$tables, "15_season_effect_means.csv"))
message("Season effect table saved.")

# ── Identify features with highest ICC_season for plotting ───────────────────

vardecomp_path <- file.path(paths$tables, "11_real_variance_decomposition.csv")
top_features <- if (file.exists(vardecomp_path)) {
  vd <- readr::read_csv(vardecomp_path, show_col_types = FALSE)
  vd |>
    dplyr::arrange(dplyr::desc(prop_b_mean)) |>
    dplyr::slice_head(n = 12) |>
    dplyr::pull(feature)
} else {
  # Fallback: features with highest posterior mean |b_j| range
  range_by_feat <- apply(B_mean, 2, function(x) diff(range(x)))
  feature_names[order(range_by_feat, decreasing = TRUE)[seq_len(12)]]
}
message("Plotting ", length(top_features), " features with highest ICC_season.")

# ── Temporal line plots (faceted grid) ───────────────────────────────────────

# Convert season labels to an ordered factor (chronological order from model)
# Use original-unit deviations for plots if available; fall back to Z-scores
use_orig <- "b_orig_mean" %in% names(season_tbl)
y_col <- if (use_orig) "b_orig_mean" else "b_mean"
ylo_col <- if (use_orig) "b_orig_q05" else "b_q05"
yhi_col <- if (use_orig) "b_orig_q95" else "b_q95"

# Determine units for y-axis label
get_unit_label <- function(feat) {
  if (startsWith(feat, "per90_")) "events per 90 min (deviation from average)"
  else if (startsWith(feat, "rate_"))  "proportion (deviation from average)"
  else                                  "value (deviation from average)"
}
y_axis_label <- if (use_orig) "Season effect in original units (deviation from average)" else "Season effect (standardised)"

plot_df <- season_tbl |>
  dplyr::filter(feature %in% top_features) |>
  dplyr::rename(y_mean = !!y_col, y_lo = !!ylo_col, y_hi = !!yhi_col) |>
  dplyr::mutate(
    season_f = factor(season, levels = season_levels),
    feature_label = gsub("per90_|rate_", "", feature),
    feature_label = gsub("_", " ", feature_label),
    feature_f = factor(feature_label,
      levels = gsub("per90_|rate_", "",
               gsub("_", " ", top_features)))
  )

groups  <- feature_group_lookup()
plot_df <- plot_df |>
  dplyr::left_join(groups, by = "feature") |>
  dplyr::mutate(group = dplyr::coalesce(group, "Other"))

p_grid <- ggplot2::ggplot(plot_df,
  ggplot2::aes(x = season_f, y = y_mean, group = feature_f, colour = group)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = y_lo, ymax = y_hi, fill = group),
    alpha = 0.18, colour = NA) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  ggplot2::scale_colour_manual(values = group_colours(), guide = "none") +
  ggplot2::scale_fill_manual(values = group_colours(), guide = "none") +
  ggplot2::facet_wrap(~feature_f, scales = "free_y", ncol = 3) +
  ggplot2::labs(
    title    = "Temporal season effects b_j — top 12 features by ICC_season",
    subtitle = "Posterior mean ± 90% credible interval",
    x        = "Season",
    y        = y_axis_label
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = ggplot2::element_text(size = 7, face = "bold")
  )
save_plot(p_grid,
  file.path(paths$figures, "15_season_profiles_grid.png"), w = 14, h = 12)

# ── Season profile heatmap (seasons × all features) ───────────────────────────

# Cluster features by pattern similarity; keep seasons in temporal order
B_plot <- B_mean
colnames(B_plot) <- feature_names
rownames(B_plot) <- season_levels

if (requireNamespace("pheatmap", quietly = TRUE)) {
  feat_ann <- data.frame(
    Group = groups$group[match(feature_names, groups$feature)],
    row.names = feature_names
  )
  ann_colors <- list(Group = group_colours()[
    intersect(names(group_colours()), unique(feat_ann$Group))])

  png(file.path(paths$figures, "15_season_profile_heatmap.png"),
      width = 1600, height = 600, res = 130)
  pheatmap::pheatmap(
    B_plot,
    cluster_rows   = FALSE,   # keep temporal order
    cluster_cols   = TRUE,
    show_rownames  = TRUE,
    show_colnames  = TRUE,
    annotation_col = feat_ann,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
    fontsize_row   = 9,
    fontsize_col   = 6,
    main           = "Posterior mean season effects b_j (all features)",
    border_color   = NA
  )
  dev.off()
  message("Saved: 15_season_profile_heatmap.png")
} else {
  heat_df <- as.data.frame(B_plot) |>
    tibble::rownames_to_column("season") |>
    tidyr::pivot_longer(-season, names_to = "feature", values_to = "effect") |>
    dplyr::left_join(groups, by = "feature") |>
    dplyr::mutate(
      season_f = factor(season, levels = season_levels),
      feature_label = gsub("per90_|rate_", "", feature),
      feature_label = gsub("_", " ", feature_label),
      group = dplyr::coalesce(group, "Other"),
      group = factor(group, levels = group_order())
    )

  p_heat <- ggplot2::ggplot(heat_df,
    ggplot2::aes(x = season_f, y = feature_label, fill = effect)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0) +
    ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      title = "Posterior mean season effects b_j",
      x = "Season", y = NULL) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, size = 7)
    )
  save_plot(p_heat,
    file.path(paths$figures, "15_season_profile_heatmap.png"), w = 14, h = 16)
}

message("15_season_profiles complete.")
