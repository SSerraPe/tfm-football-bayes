# Intraclass Correlation (ICC) analysis.
#
# The variance proportions prop_a, prop_b, prop_e computed in the Stan models
# are exactly intraclass correlations (ICCs) as derived in the professor's
# corr_diagmodel.pdf:
#
#   ICC_player[k]  = σ²_{a,k} / (σ²_{a,k} + σ²_{b,k} + σ²_{e,k})
#      = Corr(y_{ijk}, y_{ij'k})   [same player, different seasons]
#
#   ICC_season[k]  = σ²_{b,k} / (σ²_{a,k} + σ²_{b,k} + σ²_{e,k})
#      = Corr(y_{ijk}, y_{i'jk})   [different players, same season]
#
# For the low-rank model σ²_{a,k} = Σ_a[k,k] = Σ_q λ²_{kq} + ψ²_{a,k}.
# The ICC from a different player in a different season is 0 (no shared effect).
#
# This script:
#   1. Presents ICC clearly with uncertainty intervals
#   2. Compares ICC between the diagonal model (stage 03) and the low-rank
#      model (stage 10/11) to show what adding factor structure changes
#
# Reads:   outputs/tables/11_real_variance_decomposition.csv
#          outputs/tables/03_real_diagonal_additive_posterior_summary.csv
# Writes:  outputs/tables/16_icc_summary.csv
#          outputs/figures/16_icc_ranked_barchart.png
#          outputs/figures/16_icc_diagonal_vs_lowrank.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/16_icc_analysis.R"
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
P              <- length(feature_names)

save_plot <- function(p, fname, w = 12, h = 10) {
  ggplot2::ggsave(fname, p, width = w, height = h, dpi = 150)
  message("Saved: ", basename(fname))
}

# ── Load low-rank model ICC (from stage 11 variance decomposition) ────────────

lowrank_vd_path <- file.path(paths$tables, "11_real_variance_decomposition.csv")
if (!file.exists(lowrank_vd_path)) {
  message("Missing 11_real_variance_decomposition.csv — run stage 11 first.")
  quit(save = "no", status = 0)
}
lowrank_vd <- readr::read_csv(lowrank_vd_path, show_col_types = FALSE)

# Also load posterior summary for CIs on ICC
real_summary_path <- file.path(paths$tables, "10_real_lowrank_a_diag_b_posterior_summary.csv")
icc_ci <- if (file.exists(real_summary_path)) {
  ps <- readr::read_csv(real_summary_path, show_col_types = FALSE)
  pa <- extract_vector_summary(ps, "prop_a", feature_names)
  pb <- extract_vector_summary(ps, "prop_b", feature_names)
  pe <- extract_vector_summary(ps, "prop_e", feature_names)
  list(prop_a = pa, prop_b = pb, prop_e = pe)
} else NULL

# ── Build ICC summary table ───────────────────────────────────────────────────

groups <- feature_group_lookup()

icc_tbl <- lowrank_vd |>
  dplyr::rename(
    icc_player  = prop_a_mean,
    icc_season  = prop_b_mean,
    icc_residual = prop_e_mean
  ) |>
  dplyr::select(feature, icc_player, icc_season, icc_residual,
    dplyr::any_of(c("sigma_a_diag_mean", "var_b_mean", "var_e_mean"))) |>
  dplyr::left_join(groups, by = "feature") |>
  dplyr::mutate(
    group   = dplyr::coalesce(group, "Other"),
    group   = factor(group, levels = group_order())
  )

# Attach CIs from posterior summary if available
if (!is.null(icc_ci) && !is.null(icc_ci$prop_a)) {
  icc_tbl <- icc_tbl |>
    dplyr::mutate(
      icc_player_q05  = icc_ci$prop_a$q5,
      icc_player_q95  = icc_ci$prop_a$q95,
      icc_season_q05  = icc_ci$prop_b$q5,
      icc_season_q95  = icc_ci$prop_b$q95
    )
}

readr::write_csv(icc_tbl,
  file.path(paths$tables, "16_icc_summary.csv"))
message("ICC summary table saved.")

# ── Ranked ICC barchart ───────────────────────────────────────────────────────

plot_tbl <- icc_tbl |>
  dplyr::mutate(
    feature_label = gsub("per90_|rate_", "", feature),
    feature_label = gsub("_", " ", feature_label),
    feature_f     = stats::reorder(feature_label, icc_player)
  ) |>
  tidyr::pivot_longer(
    c(icc_player, icc_season, icc_residual),
    names_to = "icc_type", values_to = "icc_value"
  ) |>
  dplyr::mutate(
    icc_type = dplyr::recode(icc_type,
      icc_player   = "ICC player (same player, diff seasons)",
      icc_season   = "ICC season (diff players, same season)",
      icc_residual = "Residual"
    ),
    icc_type = factor(icc_type, levels = c(
      "ICC player (same player, diff seasons)",
      "ICC season (diff players, same season)",
      "Residual"
    ))
  )

p_bar <- ggplot2::ggplot(plot_tbl,
  ggplot2::aes(x = icc_value, y = feature_f, fill = icc_type)) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::scale_fill_manual(
    values = c(
      "ICC player (same player, diff seasons)"  = "#2166ac",
      "ICC season (diff players, same season)"  = "#f4a582",
      "Residual"                                = "#d6604d"),
    name = NULL) +
  ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_x_continuous(labels = scales::percent) +
  ggplot2::labs(
    title = "Intraclass correlations — low-rank player model (real La Liga data)",
    subtitle = paste0(
      "ICC_player = Corr(y_ijk, y_ij'k) [same player, diff seasons]\n",
      "ICC_season = Corr(y_ijk, y_i'jk) [diff players, same season]"),
    x = "Proportion of total variance", y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    strip.text.y    = ggplot2::element_text(angle = 0, hjust = 0, size = 7, face = "bold"),
    legend.position = "top",
    axis.text.y     = ggplot2::element_text(size = 7)
  )
save_plot(p_bar,
  file.path(paths$figures, "16_icc_ranked_barchart.png"), w = 11, h = 16)

# ── Comparison: diagonal model vs low-rank model ICC ─────────────────────────

diag_summary_path <- file.path(paths$tables, "03_real_diagonal_additive_posterior_summary.csv")
if (file.exists(diag_summary_path)) {
  diag_ps   <- readr::read_csv(diag_summary_path, show_col_types = FALSE)
  diag_pa   <- extract_vector_summary(diag_ps, "prop_a", feature_names)
  diag_pb   <- extract_vector_summary(diag_ps, "prop_b", feature_names)

  if (!is.null(diag_pa) && !is.null(diag_pb)) {
    compare_tbl <- tibble::tibble(
      feature          = feature_names,
      icc_player_diag  = diag_pa$mean,
      icc_season_diag  = diag_pb$mean,
      icc_player_lr    = icc_tbl$icc_player,
      icc_season_lr    = icc_tbl$icc_season
    ) |>
      dplyr::left_join(groups, by = "feature") |>
      dplyr::mutate(
        group = dplyr::coalesce(group, "Other"),
        group = factor(group, levels = group_order()),
        feature_label = gsub("per90_|rate_", "", feature),
        feature_label = gsub("_", " ", feature_label),
        feature_f     = stats::reorder(feature_label, icc_player_lr),
        delta_player  = icc_player_lr - icc_player_diag,
        delta_season  = icc_season_lr  - icc_season_diag
      )

    # Scatter: diagonal ICC_player vs low-rank ICC_player
    p_comp_player <- ggplot2::ggplot(compare_tbl,
      ggplot2::aes(x = icc_player_diag, y = icc_player_lr, colour = group)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::scale_colour_manual(values = group_colours(), name = "Group") +
      ggplot2::scale_x_continuous(labels = scales::percent) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(
        title    = "ICC_player: diagonal model vs low-rank model",
        subtitle = "Points above diagonal = low-rank model attributes more variance to players",
        x        = "ICC_player — diagonal Σ_a (stage 03)",
        y        = "ICC_player — low-rank Σ_a (stage 10)"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "right")

    p_comp_season <- ggplot2::ggplot(compare_tbl,
      ggplot2::aes(x = icc_season_diag, y = icc_season_lr, colour = group)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_point(size = 2.5, alpha = 0.85) +
      ggplot2::scale_colour_manual(values = group_colours(), name = "Group") +
      ggplot2::scale_x_continuous(labels = scales::percent) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::labs(
        title    = "ICC_season: diagonal model vs low-rank model",
        x        = "ICC_season — diagonal Σ_a (stage 03)",
        y        = "ICC_season — low-rank Σ_a (stage 10)"
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(legend.position = "right")

    # Delta plot: change in ICC_player from diagonal to low-rank, by feature
    delta_df <- compare_tbl |>
      dplyr::mutate(feature_fd = stats::reorder(feature_label, delta_player))

    p_delta <- ggplot2::ggplot(delta_df,
      ggplot2::aes(x = delta_player, y = feature_fd, fill = group)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = group_colours(), name = "Group") +
      ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_x_continuous(labels = scales::percent) +
      ggplot2::labs(
        title    = "Change in ICC_player: low-rank minus diagonal model",
        subtitle = "Positive = low-rank model assigns more variance to player effects",
        x        = "ΔICC_player", y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(
        strip.text.y    = ggplot2::element_text(angle = 0, hjust = 0, size = 7, face = "bold"),
        legend.position = "none",
        axis.text.y     = ggplot2::element_text(size = 7)
      )

    # Combined figure: scatter + delta
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p_combined <- (p_comp_player | p_comp_season) / p_delta +
        patchwork::plot_annotation(
          title = "ICC comparison: diagonal vs low-rank Σ_a model")
      save_plot(p_combined,
        file.path(paths$figures, "16_icc_diagonal_vs_lowrank.png"), w = 14, h = 14)
    } else {
      save_plot(p_comp_player,
        file.path(paths$figures, "16_icc_player_scatter.png"), w = 8, h = 6)
      save_plot(p_comp_season,
        file.path(paths$figures, "16_icc_season_scatter.png"), w = 8, h = 6)
      save_plot(p_delta,
        file.path(paths$figures, "16_icc_delta_player.png"), w = 10, h = 14)
    }
    message("ICC comparison figures saved.")
  }
} else {
  message("Stage 03 diagonal model summary not found — skipping comparison.")
}

# ── Print top/bottom ICC features ─────────────────────────────────────────────

cat("\nTop 10 features by ICC_player (most stable across seasons per player):\n")
print(icc_tbl |>
  dplyr::arrange(dplyr::desc(icc_player)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::select(feature, group, icc_player, icc_season, icc_residual))

cat("\nTop 10 features by ICC_season (most synchronised across players per season):\n")
print(icc_tbl |>
  dplyr::arrange(dplyr::desc(icc_season)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::select(feature, group, icc_player, icc_season, icc_residual))

message("16_icc_analysis complete.")
