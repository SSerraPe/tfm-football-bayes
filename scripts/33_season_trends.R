# Stage 33 — Season trend visualizations (smarter than 15_season_profiles_grid.png)
#
# Reads pre-computed 15_season_effect_means.csv (no fit re-load needed).
# Produces two complementary figures for professor_summary.qmd:
#
#   33_season_heatmap.png       — heatmap: seasons × all 48 features (overview)
#   33_season_group_trends.png  — group-averaged season effect, one line per group
#   33_season_top6.png          — line charts for top-6 features by ICC_season

source("src/bootstrap.R")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr); library(tidyr)
})
source("src/loading_visualization.R")

# ── 1. Load data ──────────────────────────────────────────────────────────────

bj   <- read_csv(file.path(paths$tables, "15_season_effect_means.csv"),
                 show_col_types = FALSE)
icc  <- read_csv(file.path(paths$tables, "16_icc_summary.csv"),
                 show_col_types = FALSE)
mo   <- readRDS(file.path(paths$processed, "model_objects.rds"))

# Wyscout season IDs → La Liga season labels (2011/12 – 2022/23)
sid_to_label <- c(
  "8504"   = "2011/12", "9630"   = "2012/13", "10837"  = "2013/14",
  "181144" = "2014/15", "185641" = "2015/16", "185828" = "2016/17",
  "186267" = "2017/18", "187526" = "2018/19", "188124" = "2019/20",
  "189012" = "2020/21", "189974" = "2021/22", "191659" = "2022/23"
)

bj <- bj |>
  mutate(
    season_lbl = sid_to_label[as.character(season)],
    season_lbl = factor(season_lbl, levels = unique(season_lbl[order(season)]))
  )

groups <- feature_group_lookup()
bj <- bj |> left_join(groups, by = "feature") |>
  mutate(group = coalesce(group, "Other"))

feature_names <- mo$variable_names
feature_label <- setNames(
  gsub("_", " ", gsub("^per90_|^rate_", "", feature_names)),
  feature_names
)

# ── 2. Heatmap: seasons × all 48 features ─────────────────────────────────────
# Group features by football category (facet rows); seasons chronological (x-axis)

heat_df <- bj |>
  mutate(
    feat_lbl = feature_label[feature],
    group    = factor(group, levels = group_order())
  ) |>
  filter(!is.na(group))

# Clamp extreme values for colour scale readability
clamp <- quantile(abs(heat_df$b_mean), 0.98, na.rm = TRUE)

p_heat <- ggplot(heat_df,
    aes(x = season_lbl, y = feat_lbl, fill = pmax(pmin(b_mean, clamp), -clamp))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "white",
    high     = "#d73027",
    midpoint = 0,
    name     = "Season\neffect",
    limits   = c(-clamp, clamp)
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Posterior mean season effects $b_{j,p}$ (stage 28, K=3)",
    subtitle = "Blue = below-average season for that feature; red = above-average",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 7),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7, face = "bold"),
    legend.position = "right",
    panel.grid   = element_blank()
  )

ggsave(file.path(paths$figures, "33_season_heatmap.png"),
       p_heat, width = 11, height = 14, dpi = 150)
message("Saved: 33_season_heatmap.png")

# ── 3. Group-average season effect ────────────────────────────────────────────

group_avg <- bj |>
  group_by(season, season_lbl, group) |>
  summarise(
    b_group_mean = mean(b_mean,  na.rm = TRUE),
    b_group_lo   = mean(b_q05,   na.rm = TRUE),
    b_group_hi   = mean(b_q95,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(group = factor(group, levels = group_order()))

p_group <- ggplot(group_avg,
    aes(x = season_lbl, y = b_group_mean, colour = group, group = group)) +
  geom_ribbon(aes(ymin = b_group_lo, ymax = b_group_hi, fill = group),
              alpha = 0.12, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = group_colours(), name = "Feature group") +
  scale_fill_manual(  values = group_colours(), name = "Feature group") +
  labs(
    title    = "Group-averaged season effects $\\hat{b}_{j}$ (stage 28, K=3)",
    subtitle = "Each line = mean posterior season effect averaged within a feature group. Ribbon = mean of 90% CIs.",
    x = "Season", y = "Mean season effect (std. units)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(paths$figures, "33_season_group_trends.png"),
       p_group, width = 11, height = 6, dpi = 150)
message("Saved: 33_season_group_trends.png")

# ── 4. Top-6 individual features by ICC_season ────────────────────────────────

top6 <- icc |>
  filter(!is.na(icc_season)) |>
  arrange(desc(icc_season)) |>
  slice_head(n = 6) |>
  pull(feature)

cat("Top-6 features by ICC_season:\n")
print(top6)

top6_df <- bj |>
  filter(feature %in% top6) |>
  mutate(
    feat_lbl = feature_label[feature],
    icc_s    = icc$icc_season[match(feature, icc$feature)],
    feat_f   = factor(feat_lbl,
      levels = feature_label[top6[order(icc$icc_season[match(top6, icc$feature)],
                                        decreasing = TRUE)]])
  )

p_top6 <- ggplot(top6_df,
    aes(x = season_lbl, y = b_mean, colour = group, group = feat_f)) +
  geom_ribbon(aes(ymin = b_q05, ymax = b_q95, fill = group),
              alpha = 0.18, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = group_colours(), guide = "none") +
  scale_fill_manual(  values = group_colours(), guide = "none") +
  facet_wrap(~feat_f, scales = "free_y", ncol = 2) +
  labs(
    title    = "Top-6 features by ICC_season: posterior season effects $\\hat{b}_{j,p}$",
    subtitle = "Ribbon = 90% credible interval. Features ordered by ICC_season (decreasing).",
    x = "Season", y = "Season effect (std. units)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text  = element_text(size = 9, face = "bold")
  )

ggsave(file.path(paths$figures, "33_season_top6.png"),
       p_top6, width = 10, height = 9, dpi = 150)
message("Saved: 33_season_top6.png")

message("\nStage 33 complete.")
