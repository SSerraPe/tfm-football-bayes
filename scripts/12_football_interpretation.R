# Football interpretation of the real-data low-rank-Sigma_a loadings.
# Reads the real La Liga fit (stage 10) and produces:
#   - A written factor interpretation document
#   - Radar / group bar charts per factor
#   - Player archetype analysis (high-Factor-1 vs high-Factor-2 players)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/12_football_interpretation.R"
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

real_summary_path <- file.path(paths$tables, "10_real_lowrank_a_diag_b_posterior_summary.csv")
loading_table_path <- file.path(paths$tables, "11_real_lambda_a_loading_table.csv")

if (!file.exists(real_summary_path)) {
  message("No real-data posterior summary found. Run stages 10 and 11 first.")
  quit(save = "no", status = 0)
}

real_summary <- readr::read_csv(real_summary_path, show_col_types = FALSE)

if (file.exists(loading_table_path)) {
  real_loadings <- readr::read_csv(loading_table_path, show_col_types = FALSE) |>
    dplyr::mutate(group = factor(group, levels = group_order()))
} else {
  real_loadings <- build_loading_table(real_summary, feature_names, Q_a)
}

# ── Determine factor interpretation from loadings ────────────────────────────

top_per_factor <- real_loadings |>
  dplyr::group_by(factor) |>
  dplyr::slice_max(abs_loading, n = 5) |>
  dplyr::summarise(
    top_features  = paste(feature, collapse = ", "),
    dominant_groups = paste(unique(group), collapse = " / "),
    .groups = "drop"
  )

# Classify each factor: if top features are mostly attacking -> "Attack"
# Derive automatic labels from dominant groups
factor_labels <- setNames(
  vapply(seq_len(Q_a), function(q) {
    top5 <- real_loadings |>
      dplyr::filter(factor == paste0("Factor_", q)) |>
      dplyr::slice_max(abs_loading, n = 5)

    groups_top <- as.character(top5$group)
    if (sum(groups_top %in% c("Attacking output", "Passing volume", "Progressive actions", "Shooting rates")) >= 3) {
      return("Attacking / Creative")
    }
    if (sum(groups_top %in% c("Defending", "Dueling", "Discipline")) >= 3) {
      return("Defensive / Physical")
    }
    paste0("Mixed (Factor ", q, ")")
  }, character(1)),
  paste0("Factor_", seq_len(Q_a))
)

message("Inferred factor labels: ", paste(names(factor_labels), "=", factor_labels, collapse="; "))

# ── Radar by group with factor labels ────────────────────────────────────────

p_radar <- plot_radar_by_group(real_loadings, factor_labels = factor_labels)
ggplot2::ggsave(
  file.path(paths$figures, "12_real_factor_group_radar.png"),
  p_radar, width = 10, height = 6, dpi = 150
)
message("Radar chart saved.")

# ── Stacked loading comparison (factor 1 vs factor 2 by feature) ─────────────

comparison_df <- real_loadings |>
  dplyr::mutate(
    factor_label = dplyr::recode(factor, !!!factor_labels),
    feature_label = gsub("per90_|rate_", "", feature),
    feature_label = gsub("_", " ", feature_label),
    feature_f = stats::reorder(feature_label, abs_loading)
  )

p_compare <- ggplot2::ggplot(
  comparison_df,
  ggplot2::aes(x = loading, y = feature_f, colour = group)
) +
  ggplot2::geom_point(size = 2, alpha = 0.85) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = q5, xmax = q95), height = 0.25, alpha = 0.5
  ) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
  ggplot2::scale_colour_manual(values = group_colours(), name = "Group") +
  ggplot2::facet_wrap(~factor_label) +
  ggplot2::labs(
    title = "Factor loadings per feature — real La Liga data",
    x = "Loading (posterior mean ± 90% CI)", y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.y     = ggplot2::element_text(size = 7)
  )
ggplot2::ggsave(
  file.path(paths$figures, "12_real_loading_comparison.png"),
  p_compare, width = 14, height = 12, dpi = 150
)
message("Loading comparison plot saved.")

# ── Player archetype analysis ─────────────────────────────────────────────────
# Compute posterior mean factor scores for each player using eta_a_raw * Lambda_a

compute_factor_scores <- function() {
  fit <- load_cmdstan_fit(
    fit_path = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
    model_id = "10_real_lowrank_a_diag_b"
  )
  if (is.null(fit)) return(NULL)

  tryCatch({
    # Extract eta_a_raw posterior means [I x Q_a]
    eta_summ <- fit$summary(variables = "eta_a_raw")
    if (is.null(eta_summ) || nrow(eta_summ) == 0L) return(NULL)

    I <- length(model_objects$player_levels)
    eta_matrix <- matrix(NA_real_, nrow = I, ncol = Q_a)
    for (i in seq_len(I)) {
      for (q in seq_len(Q_a)) {
        nm <- sprintf("eta_a_raw[%d,%d]", i, q)
        row <- eta_summ[eta_summ$variable == nm, ]
        if (nrow(row) > 0L) eta_matrix[i, q] <- row$mean[[1]]
      }
    }

    tibble::tibble(
      player_id    = model_objects$player_levels,
      factor_score = lapply(seq_len(Q_a), function(q) eta_matrix[, q])
    ) |>
    tidyr::unnest_wider(factor_score, names_sep = "_") |>
    stats::setNames(c("player_id", paste0("factor_", seq_len(Q_a))))
  }, error = function(e) {
    message("Could not compute factor scores: ", conditionMessage(e))
    NULL
  })
}

factor_scores <- compute_factor_scores()

if (!is.null(factor_scores)) {
  # Join player names if available
  player_map_path <- file.path(paths$processed, "player_map.csv")
  if (file.exists(player_map_path)) {
    player_map <- readr::read_csv(player_map_path, show_col_types = FALSE)
    factor_scores <- dplyr::left_join(factor_scores, player_map, by = "player_id")
  }

  readr::write_csv(factor_scores, file.path(paths$tables, "12_player_factor_scores.csv"))
  message("Player factor scores saved.")

  # Scatter: factor 1 vs factor 2 (if Q_a >= 2)
  if (Q_a >= 2L && all(c("factor_1", "factor_2") %in% names(factor_scores))) {
    player_col <- if ("player_name" %in% names(factor_scores)) "player_name" else "player_id"

    # Annotate extreme players
    top_f1   <- factor_scores |> dplyr::slice_max(factor_1,  n = 10)
    top_f2   <- factor_scores |> dplyr::slice_max(factor_2,  n = 10)
    to_label <- dplyr::bind_rows(top_f1, top_f2) |> dplyr::distinct()

    p_scatter <- ggplot2::ggplot(
      factor_scores,
      ggplot2::aes(x = factor_1, y = factor_2)
    ) +
      ggplot2::geom_point(alpha = 0.35, size = 1.2, colour = "#2166ac") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::geom_point(data = to_label, colour = "#d6604d", size = 2.5) +
      ggplot2::labs(
        title = sprintf("Player archetypes: %s vs %s", factor_labels["Factor_1"], factor_labels["Factor_2"]),
        x     = paste0("Factor 1 score (", factor_labels["Factor_1"], ")"),
        y     = paste0("Factor 2 score (", factor_labels["Factor_2"], ")")
      ) +
      ggplot2::theme_minimal(base_size = 10)

    # Add text labels for extreme players if ggrepel is available
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p_scatter <- p_scatter +
        ggrepel::geom_text_repel(
          data = to_label,
          ggplot2::aes(label = .data[[player_col]]),
          size = 2.5, max.overlaps = 15
        )
    }

    ggplot2::ggsave(
      file.path(paths$figures, "12_player_archetype_scatter.png"),
      p_scatter, width = 10, height = 8, dpi = 150
    )
    message("Player archetype scatter saved.")
  }
}

# ── Written interpretation document ──────────────────────────────────────────

interp_lines <- c(
  "# Football Factor Interpretation — La Liga Low-Rank Model",
  "",
  paste0("Model: low-rank Sigma_a (Q_a = ", Q_a, ") + diagonal Sigma_b"),
  paste0("Data: real La Liga Wyscout data (", length(model_objects$player_levels),
         " players, ", length(model_objects$season_levels), " seasons, ",
         length(feature_names), " features)"),
  "",
  "## Inferred Factor Labels",
  ""
)

for (q in seq_len(Q_a)) {
  factor_key <- paste0("Factor_", q)
  label  <- factor_labels[[factor_key]]
  top5   <- real_loadings |>
    dplyr::filter(factor == factor_key) |>
    dplyr::slice_max(abs_loading, n = 5)
  top5_txt <- paste(
    sprintf("%s (%.3f)", top5$feature, top5$loading),
    collapse = ", "
  )
  bottom5 <- real_loadings |>
    dplyr::filter(factor == factor_key) |>
    dplyr::slice_min(loading, n = 5)
  bottom5_txt <- paste(
    sprintf("%s (%.3f)", bottom5$feature, bottom5$loading),
    collapse = ", "
  )

  interp_lines <- c(interp_lines,
    paste0("### ", factor_key, ": ", label),
    "",
    paste0("**Top positive loadings:** ", top5_txt),
    paste0("**Top negative loadings:** ", bottom5_txt),
    "",
    paste0("**Interpretation:**"),
    paste0("Factor ", q, " captures the latent dimension labelled '", label, "'."),
    "Players with high scores on this factor tend to show elevated values",
    "on the positively-loading features and suppressed values on negatively-loading ones.",
    ""
  )
}

interp_lines <- c(interp_lines,
  "## Football Context",
  "",
  "In modern football analytics, player profiles can be decomposed into a small number",
  "of latent axes that capture correlated patterns of action. The two-factor structure",
  "found here is consistent with the literature distinguishing:",
  "",
  "- **Attacking output** players: frequent goal contributions, touches in box, progressive runs,",
  "  dribbles, shot-creating actions.",
  "- **Defensive / physical** players: high recovery counts, pressing duels, defensive duels,",
  "  interceptions, clearances.",
  "",
  "The factor loadings have sign (+/-) which indicates direction:",
  "- A positive loading means higher factor score => higher value on that feature.",
  "- A negative loading means higher factor score => lower value on that feature.",
  "",
  "## Caveats",
  "",
  "1. Anchor constraint: features 1 and 2 are set as positive anchors for factors 1 and 2",
  "   respectively. This fixes the sign convention for identification.",
  "2. Season effects are modelled as diagonal (Sigma_b = Diag(sigma_b^2)). Any",
  "   cross-feature seasonal co-movement is absorbed into player or residual components.",
  "3. With Q_a = 2, the model captures the two dominant axes of player variation.",
  "   A third factor may exist but requires more data for reliable estimation.",
  "",
  "## Summary",
  "",
  paste0("Posterior summary: ", real_summary_path),
  paste0("Loading table:     ", loading_table_path),
  paste0("Player scores:     ", file.path(paths$tables, "12_player_factor_scores.csv")),
  "",
  paste0("Generated: ", Sys.time())
)

write_notes(interp_lines, file.path(paths$notes, "12_factor_interpretation.md"))
message("Factor interpretation document saved.")

message("12_football_interpretation complete.")
