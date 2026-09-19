# Stage 45 — Top-5/top-5 ICC stacked-bar figure for the Results chapter (thesis
# revision Phase 3, replacing the 15-row "highest player ICC" table in
# results.tex Sec. 4 with a figure that carries the player-dominance /
# season-drift message visually rather than through 15 rows of CIs).
#
# Two small panels: the 5 features with the highest posterior-median ICC_player,
# and the 5 features with the highest posterior-median ICC_season -- each bar a
# player/season/residual stacked share summing to 100%. No overlap between the
# two feature sets in the current K=3 fit.
#
# Reads:  outputs/tables/41_icc_summary_k3.csv (already computed, stage 41)
# Writes: outputs/figures/45_icc_top5_stacked.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/45_icc_top5_barplot.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

icc <- read_csv(file.path(paths$tables, "41_icc_summary_k3.csv"), show_col_types = FALSE)

clean_name <- function(x) tools::toTitleCase(gsub("_", " ", gsub("^per90_|^rate_", "", x)))

top_player <- icc |> arrange(desc(icc_player)) |> slice_head(n = 5) |> mutate(panel = "Highest player ICC")
top_season <- icc |> arrange(desc(icc_season)) |> slice_head(n = 5) |> mutate(panel = "Highest season ICC")

plot_df <- bind_rows(top_player, top_season) |>
  mutate(feature_label = clean_name(feature)) |>
  select(panel, feature_label, icc_player, icc_season, icc_residual) |>
  pivot_longer(cols = c(icc_player, icc_season, icc_residual), names_to = "component", values_to = "share") |>
  mutate(
    component = factor(component,
      levels = c("icc_residual", "icc_season", "icc_player"),
      labels = c("Residual", "Season", "Player")
    ),
    panel = factor(panel, levels = c("Highest player ICC", "Highest season ICC"))
  )

# Order bars within each panel by that panel's own sort key (descending), independently.
order_df <- bind_rows(
  top_player |> transmute(panel = "Highest player ICC", feature_label = clean_name(feature), sort_key = icc_player),
  top_season |> transmute(panel = "Highest season ICC", feature_label = clean_name(feature), sort_key = icc_season)
)
plot_df <- plot_df |>
  left_join(order_df, by = c("panel", "feature_label")) |>
  mutate(feature_label = reorder(feature_label, sort_key))

pal <- c("Player" = "#2b81ad", "Season" = "#e08a2e", "Residual" = "#c9c9c9")

p <- ggplot(plot_df, aes(x = feature_label, y = share, fill = component)) +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~panel, scales = "free_y") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Share of total variance") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(paths$figures, "45_icc_top5_stacked.png"), p, width = 9, height = 3.6, dpi = 150)
message("Stage 45 complete: outputs/figures/45_icc_top5_stacked.png")
