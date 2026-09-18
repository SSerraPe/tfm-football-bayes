# Stage 44 — Posterior-uncertainty scatter: player positions are not points.
#
# Every player-similarity figure so far (archetype scatter, style-vs-level scatter)
# plots a single point per player (the posterior mean). This makes the point
# results.tex already argues in prose visible directly: a player's position in
# style space has a posterior *cloud* around it, and whether two players are
# "confidently close" or "ambiguously close" is a visual fact about how much
# their clouds overlap -- not just a single distance number.
#
# Two panels, both illustrating findings already written up in results.tex §8.3:
#   A) Messi + Neymar: the confident pair (closest neighbour stable in 92% of draws)
#   B) Ramos + Piqué + Umtiti + Javi García: the near four-way statistical tie
#      (each the single closest neighbour in only 16-24% of draws)
#
# Reads:  fits/38_style_score_draws.rds  (cached: draws x I x K, from stage 38 --
#         no new CmdStan extraction needed)
#         outputs/tables/29_player_factor_scores.csv  (player_index lookup)
# Writes: outputs/figures/44_posterior_clouds.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/44_posterior_clouds.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

F_draws <- readRDS(file.path(paths$fits, "38_style_score_draws.rds"))  # n_draws x I x K
n_draws <- dim(F_draws)[1]

players <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"), show_col_types = FALSE)

# Standardise PC1/PC2 across the population, per draw, matching stage 38/42's convention,
# then subsample draws for a legible scatter (4000 points/player is too dense to read).
set.seed(20260918L)
plot_idx <- sample(seq_len(n_draws), 600L)

cloud_for <- function(player_index, player_name, group) {
  pts <- lapply(plot_idx, function(s) {
    Fs_std <- scale(F_draws[s, , 1:2])
    Fs_std[player_index, ]
  })
  m <- do.call(rbind, pts)
  tibble(PC1 = m[, 1], PC2 = m[, 2], player_name = player_name, group = group)
}

panel_a <- bind_rows(
  cloud_for(63,  "L. Messi", "A: confident pair"),
  cloud_for(583, "Neymar",   "A: confident pair")
)
panel_b <- bind_rows(
  cloud_for(33,  "Sergio Ramos", "B: near four-way tie"),
  cloud_for(51,  "Gerard Piqué", "B: near four-way tie"),
  cloud_for(512, "S. Umtiti",    "B: near four-way tie"),
  cloud_for(613, "Javi García",  "B: near four-way tie")
)
plot_df <- bind_rows(panel_a, panel_b)

panel_labels <- c(
  "A: confident pair" = "A: a confident pair (Messi-Neymar, closest in 92% of draws)",
  "B: near four-way tie" = "B: a near four-way tie (each closest in only 16-24% of draws)"
)
plot_df$panel_label <- panel_labels[plot_df$group]

p <- ggplot(plot_df, aes(PC1, PC2, color = player_name)) +
  geom_point(size = 0.6, alpha = 0.25) +
  stat_ellipse(aes(group = player_name), level = 0.9, linewidth = 0.9, alpha = 0.9) +
  facet_wrap(~panel_label, scales = "free") +
  labs(
    title = "Posterior clouds, not points: two players' style position, drawn from 600 posterior draws each",
    subtitle = "Solid line = 90% posterior ellipse. Heavy overlap = a genuinely ambiguous 'closest neighbour' claim.",
    x = "PC1 (standardised, per draw)", y = "PC2 (standardised, per draw)", color = "Player"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 9))

ggsave(file.path(paths$figures, "44_posterior_clouds.png"), p, width = 11, height = 5.5, dpi = 150)
message("Saved: 44_posterior_clouds.png")
