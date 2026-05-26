# Helpers for loading table display and visualization.

feature_group_lookup <- function() {
  tibble::tribble(
    ~feature,                                   ~group,
    "per90_goals",                              "Attacking output",
    "per90_shots",                              "Attacking output",
    "per90_xg_shot",                            "Attacking output",
    "per90_head_shots",                         "Attacking output",
    "per90_touch_in_box",                       "Attacking output",
    "per90_assists",                            "Attacking output",
    "per90_xg_assist",                          "Attacking output",
    "rate_shots_on_target",                     "Shooting rates",
    "rate_goal_conversion",                     "Shooting rates",
    "per90_received_pass",                      "Passing volume",
    "per90_passes",                             "Passing volume",
    "per90_forward_passes",                     "Passing volume",
    "per90_back_passes",                        "Passing volume",
    "per90_lateral_passes",                     "Passing volume",
    "per90_long_passes",                        "Passing volume",
    "per90_smart_passes",                       "Passing volume",
    "per90_progressive_passes",                 "Passing volume",
    "per90_passes_to_final_third",              "Passing volume",
    "per90_key_passes",                         "Passing volume",
    "per90_through_passes",                     "Passing volume",
    "per90_crosses",                            "Passing volume",
    "per90_shot_assists",                       "Passing volume",
    "rate_successful_passes",                   "Passing quality",
    "rate_successful_progressive_passes",       "Passing quality",
    "rate_successful_crosses",                  "Passing quality",
    "per90_dribbles",                           "Progressive actions",
    "per90_progressive_run",                    "Progressive actions",
    "per90_accelerations",                      "Progressive actions",
    "rate_successful_dribbles",                 "Progressive actions",
    "per90_losses",                             "Defending",
    "per90_own_half_losses",                    "Defending",
    "per90_dangerous_own_half_losses",          "Defending",
    "per90_recoveries",                         "Defending",
    "per90_opponent_half_recoveries",           "Defending",
    "per90_dangerous_opponent_half_recoveries", "Defending",
    "per90_interceptions",                      "Defending",
    "per90_clearances",                         "Defending",
    "per90_shots_blocked",                      "Defending",
    "per90_sliding_tackles",                    "Defending",
    "per90_dribbles_against",                   "Defending",
    "rate_defensive_duels_won",                 "Defending",
    "rate_dribbles_against_won",                "Defending",
    "rate_successful_sliding_tackles",          "Defending",
    "per90_defensive_duels",                    "Dueling",
    "per90_offensive_duels",                    "Dueling",
    "per90_aerial_duels",                       "Dueling",
    "per90_pressing_duels",                     "Dueling",
    "per90_loose_ball_duels",                   "Dueling",
    "rate_offensive_duels_won",                 "Dueling",
    "rate_aerial_duels_won",                    "Dueling",
    "per90_fouls",                              "Discipline",
    "per90_fouls_suffered",                     "Discipline",
    "per90_yellow_cards",                       "Discipline"
  )
}

group_order <- function() {
  c("Attacking output", "Shooting rates", "Passing volume", "Passing quality",
    "Progressive actions", "Defending", "Dueling", "Discipline")
}

group_colours <- function() {
  c(
    "Attacking output"   = "#E63946",
    "Shooting rates"     = "#FF6B6B",
    "Passing volume"     = "#2196F3",
    "Passing quality"    = "#64B5F6",
    "Progressive actions"= "#9C27B0",
    "Defending"          = "#4CAF50",
    "Dueling"            = "#FF9800",
    "Discipline"         = "#795548"
  )
}

build_loading_table <- function(posterior_summary, feature_names, Q_a) {
  check_packages(c("dplyr", "tibble", "stringr"))

  groups <- feature_group_lookup()

  rows <- lapply(seq_len(Q_a), function(q) {
    means <- vapply(seq_along(feature_names), function(p) {
      nm <- sprintf("Lambda_a[%d,%d]", p, q)
      row <- posterior_summary[posterior_summary$variable == nm, ]
      if (nrow(row) == 0L) return(NA_real_)
      row$mean[[1]]
    }, numeric(1))

    q5 <- vapply(seq_along(feature_names), function(p) {
      nm <- sprintf("Lambda_a[%d,%d]", p, q)
      row <- posterior_summary[posterior_summary$variable == nm, ]
      if (nrow(row) == 0L) return(NA_real_)
      row$q5[[1]]
    }, numeric(1))

    q95 <- vapply(seq_along(feature_names), function(p) {
      nm <- sprintf("Lambda_a[%d,%d]", p, q)
      row <- posterior_summary[posterior_summary$variable == nm, ]
      if (nrow(row) == 0L) return(NA_real_)
      row$q95[[1]]
    }, numeric(1))

    tibble::tibble(
      feature  = feature_names,
      factor   = paste0("Factor_", q),
      loading  = means,
      q5       = q5,
      q95      = q95
    )
  })

  tbl <- dplyr::bind_rows(rows) |>
    dplyr::left_join(groups, by = "feature") |>
    dplyr::mutate(
      group = dplyr::coalesce(group, "Other"),
      group = factor(group, levels = group_order()),
      sign  = dplyr::if_else(loading >= 0, "+", "-"),
      abs_loading = abs(loading),
      communality = NA_real_
    )

  # Add communality (sum of squared loadings across factors per feature)
  comm <- tbl |>
    dplyr::group_by(feature) |>
    dplyr::summarise(communality = sum(loading^2, na.rm = TRUE), .groups = "drop")
  tbl <- tbl |>
    dplyr::select(-communality) |>
    dplyr::left_join(comm, by = "feature")

  tbl
}

plot_loading_heatmap <- function(loading_table, title = "Factor Loadings (Λ_a)") {
  check_packages(c("ggplot2", "dplyr"))

  df <- loading_table |>
    dplyr::arrange(group, feature) |>
    dplyr::mutate(
      feature_label = gsub("per90_|rate_", "", feature),
      feature_label = gsub("_", " ", feature_label),
      feature_f = factor(feature_label, levels = rev(unique(feature_label)))
    )

  ggplot2::ggplot(df, ggplot2::aes(x = factor, y = feature_f, fill = loading)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", loading)),
      size = 2.5, colour = "white", fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d6604d",
      midpoint = 0, limits = c(-1, 1), oob = scales::squish,
      name = "Loading"
    ) +
    ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      title = title,
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, size = 8, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 7),
      panel.grid   = ggplot2::element_blank(),
      legend.position = "right"
    )
}

plot_top_loadings <- function(loading_table, top_n = 12) {
  check_packages(c("ggplot2", "dplyr"))

  df <- loading_table |>
    dplyr::group_by(factor) |>
    dplyr::slice_max(abs_loading, n = top_n) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      feature_label = gsub("per90_|rate_", "", feature),
      feature_label = gsub("_", " ", feature_label)
    )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = loading,
      y = stats::reorder(paste0(feature_label, " (", factor, ")"), loading),
      fill = group
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = q5, xmax = q95),
      height = 0.3, colour = "grey40"
    ) +
    ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = group_colours(), name = "Group") +
    ggplot2::facet_wrap(~factor, scales = "free_y") +
    ggplot2::labs(
      title = "Top loadings by factor (posterior mean ± 90% CI)",
      x = "Loading", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
}

plot_radar_by_group <- function(loading_table, factor_labels = NULL) {
  check_packages(c("ggplot2", "dplyr", "tidyr"))

  avg_by_group <- loading_table |>
    dplyr::group_by(factor, group) |>
    dplyr::summarise(mean_abs_loading = mean(abs_loading, na.rm = TRUE), .groups = "drop")

  if (!is.null(factor_labels)) {
    avg_by_group <- avg_by_group |>
      dplyr::mutate(factor_label = dplyr::recode(factor, !!!factor_labels))
  } else {
    avg_by_group <- dplyr::mutate(avg_by_group, factor_label = factor)
  }

  ggplot2::ggplot(
    avg_by_group,
    ggplot2::aes(x = group, y = mean_abs_loading, fill = factor_label, group = factor_label)
  ) +
    ggplot2::geom_col(position = "dodge", alpha = 0.85) +
    ggplot2::scale_fill_brewer(palette = "Set1", name = "Factor") +
    ggplot2::labs(
      title = "Mean |loading| by feature group and factor",
      x = NULL, y = "Mean |loading|"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "top"
    )
}
