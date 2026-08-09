# Stage 36 — Top 5 most-changed season profiles
# Finds the 5 features with the largest peak-to-trough swing in season effects
# and produces one panel per feature (line + 90% CI ribbon).

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "scripts/36_season_top5_profiles.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})
source(file.path(model_root, "src", "loading_visualization.R"))

# ── 1. Load data ──────────────────────────────────────────────────────────────
season_eff <- read_csv(file.path(paths$tables, "15_season_effect_means.csv"),
                       show_col_types = FALSE)
season_map <- read_csv(file.path(paths$processed, "season_map.csv"),
                       show_col_types = FALSE)

# Season ID → readable label (e.g. "11/12")
# season_map: season_index (1..12), season_id (8504, 9630, …)
season_label_vec <- c(
  "8504"   = "11/12", "9630"   = "12/13", "10837"  = "13/14",
  "181144" = "14/15", "185641" = "15/16", "185828" = "16/17",
  "186267" = "17/18", "187526" = "18/19", "188124" = "19/20",
  "189012" = "20/21", "189974" = "21/22", "191659" = "22/23"
)

season_eff <- season_eff |>
  mutate(season_label = season_label_vec[as.character(season)],
         season_label = factor(season_label, levels = season_label_vec))

# ── 2. Pick top-5 features by range of b_mean ─────────────────────────────────
feature_ranges <- season_eff |>
  group_by(feature) |>
  summarise(
    b_range = max(b_mean) - min(b_mean),
    b_max   = max(b_mean),
    b_min   = min(b_mean),
    .groups = "drop"
  ) |>
  arrange(desc(b_range))

top5_features <- feature_ranges$feature[1:5]

cat("Top 5 features by season-effect range (max b_mean - min b_mean):\n")
print(as.data.frame(feature_ranges[1:5, c("feature", "b_range", "b_max", "b_min")]))

# ── 3. Group colour lookup ────────────────────────────────────────────────────
groups <- feature_group_lookup()
gcols  <- group_colours()

# ── 4. Build one panel per top-5 feature ─────────────────────────────────────
make_panel <- function(feat, rank_i) {
  df   <- season_eff |> filter(feature == feat)
  grp  <- groups$group[groups$feature == feat]
  if (length(grp) == 0) grp <- "Other"
  col  <- if (grp %in% names(gcols)) gcols[[grp]] else "#555555"

  feat_label <- gsub("_", " ", gsub("^per90_|^rate_", "", feat))
  feat_label <- paste0(toupper(substr(feat_label, 1, 1)),
                       substr(feat_label, 2, nchar(feat_label)))
  rng_val    <- round(max(df$b_mean) - min(df$b_mean), 3)
  peak_season <- df$season_label[which.max(df$b_mean)]
  trough_season <- df$season_label[which.min(df$b_mean)]

  ggplot(df, aes(x = season_label, y = b_mean, group = 1)) +
    geom_ribbon(aes(ymin = b_q05, ymax = b_q95),
                fill = col, alpha = 0.20, colour = NA) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_line(colour = col, linewidth = 1.1) +
    geom_point(colour = col, size = 2.2) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("Δ = %.3f σ\npeak: %s | trough: %s",
                             rng_val, peak_season, trough_season),
             hjust = 1.05, vjust = 1.4, size = 2.9, colour = "grey40",
             lineheight = 1.3) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title    = sprintf("#%d  %s", rank_i, feat_label),
      subtitle = grp,
      x        = NULL,
      y        = "Season effect (σ)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 10),
      plot.subtitle    = element_text(size = 8, colour = "grey40"),
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.minor = element_blank()
    )
}

panels <- mapply(make_panel, top5_features, seq_along(top5_features),
                 SIMPLIFY = FALSE)

# ── 5. Combine with patchwork ─────────────────────────────────────────────────
fig_combined <- wrap_plots(panels, ncol = 1) +
  plot_annotation(
    title    = "Top 5 features by season-effect range",
    subtitle = paste(
      "Posterior mean season effect b̂ⱼ with 90% credible intervals.",
      "Range = max − min of posterior mean across 12 La Liga seasons (2011/12–2022/23)."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

out_fig <- file.path(paths$figures, "36_season_top5_profiles.png")
ggsave(out_fig, fig_combined, width = 9, height = 14, dpi = 150)
cat("Saved:", out_fig, "\n")

# ── 6. Text summary ───────────────────────────────────────────────────────────
direction_label <- function(feat) {
  df   <- season_eff |> filter(feature == feat) |> arrange(season_label)
  vals <- df$b_mean
  if (all(diff(vals) >= -1e-9)) return("rising")
  if (all(diff(vals) <=  1e-9)) return("falling")
  return("non-monotone")
}

summary_lines <- c(
  "=== Stage 36: Top 5 features by season-effect range ===",
  sprintf("Date: %s", Sys.Date()),
  "",
  "Measure: max(b_mean) - min(b_mean) across 12 La Liga seasons.",
  "",
  sprintf("%-55s %9s %13s", "Feature", "Range (σ)", "Direction")
)
for (feat in top5_features) {
  df_f <- season_eff |> filter(feature == feat)
  rng  <- round(max(df_f$b_mean) - min(df_f$b_mean), 3)
  dir  <- direction_label(feat)
  summary_lines <- c(summary_lines,
                     sprintf("%-55s %9.3f %13s", feat, rng, dir))
}

out_note <- file.path(paths$notes, "36_season_top5_summary.txt")
writeLines(summary_lines, out_note)
cat("Saved:", out_note, "\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
