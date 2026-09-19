# Stage 50 -- Forest plot of canonical factor scores under the fixed reference
# orientation (Results-chapter rewrite).
#
# 3x3 grid: rows = the three case-study targets (Messi, Ramos, Xavi), columns = PC1/PC2/
# PC3. Each cell holds 7 horizontal 90% credible intervals, same order in every cell:
# target, its 3 style-neighbours, its 3 most distant players (rule and picks documented in
# stage 48). Built as three manually-positioned (non-faceted) column plots combined with
# patchwork, rather than facet_wrap/facet_grid, so that the target's own name is the only
# thing marking each row block (no redundant rotated side label or repeated strip text) and
# the x-scale can be identical, symmetric, and shared globally without free-scale faceting
# workarounds.
#
# Reads:  outputs/tables/49_score_credible_intervals.csv (player x pc x {q05,q50,q95}
#           under the fixed reference orientation)
#         outputs/tables/48_distant_players.csv, outputs/tables/37_player_similarity_top10.csv
#           (row order/composition per target)
# Writes: images/50_score_forest_plot.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/50_score_forest_plot.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "patchwork"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

ci      <- read_csv(file.path(paths$tables, "49_score_credible_intervals.csv"), show_col_types = FALSE)
distant <- read_csv(file.path(paths$tables, "48_distant_players.csv"), show_col_types = FALSE)
top10   <- read_csv(file.path(paths$tables, "37_player_similarity_top10.csv"), show_col_types = FALSE)

TARGETS <- c("L. Messi", "Sergio Ramos", "Xavi")
TARGET_COL <- c("L. Messi" = "#C0392B", "Sergio Ramos" = "#1F6F8B", "Xavi" = "#7B3294")
DISTANT_COL <- "#7A7A7A"

# ---- build the 7-row-per-target row order, with a "role" tag ---------------------------
row_order <- lapply(TARGETS, function(tg) {
  nb <- top10 |> filter(player_name == tg, rank <= 3) |> arrange(rank) |> pull(similar_player)
  di <- distant |> filter(target == tg) |> arrange(rank) |> pull(player_name)
  tibble(target = tg, player_name = c(tg, nb, di),
         role = c("target", rep("neighbour", 3), rep("distant", 3)),
         row_in_block = 1:7)
}) |> bind_rows()

# y positions: 7 rows per target block, one blank row of gap between blocks, target 1's
# block at the top.
block_gap <- 1
row_order <- row_order |>
  group_by(target) |>
  mutate(block_idx = match(target, TARGETS)) |>
  ungroup() |>
  mutate(y = -( (block_idx - 1) * (7 + block_gap) + row_in_block ))

plot_df <- row_order |>
  left_join(ci, by = "player_name", relationship = "many-to-many") |>
  mutate(
    colour = ifelse(role == "distant", DISTANT_COL, TARGET_COL[target]),
    shape  = case_when(role == "target" ~ 23, role == "neighbour" ~ 21, role == "distant" ~ 22),
    size   = case_when(role == "target" ~ 4.2, role == "neighbour" ~ 3.0, role == "distant" ~ 3.0)
  )

# Global symmetric x-limits, shared across every panel (a superset of "shared within a
# column" -- simpler, and the three axes are on directly comparable canonical-score scales).
xmax <- max(abs(c(plot_df$q05, plot_df$q95))) * 1.08

make_col <- function(pcx, show_y_labels) {
  df <- plot_df |> filter(pc == pcx)
  labels_df <- row_order |> filter(!duplicated(paste(target, y)))
  p <- ggplot(df, aes(y = y)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "#AAAAAA") +
    geom_segment(aes(x = q05, xend = q95, yend = y, colour = colour),
                 linewidth = 1.5, alpha = 0.85, lineend = "round") +
    geom_point(aes(x = q50, fill = colour, shape = shape, size = size),
               colour = "white", stroke = 0.7) +
    scale_colour_identity() + scale_fill_identity() + scale_shape_identity() + scale_size_identity() +
    scale_x_continuous(limits = c(-xmax, xmax)) +
    scale_y_continuous(breaks = labels_df$y,
                        labels = if (show_y_labels) labels_df$player_name else NULL,
                        expand = expansion(add = 0.6)) +
    labs(x = NULL, y = NULL, title = pcx) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(family = "serif", size = 7.8, hjust = 1),
      axis.text.x = element_text(family = "serif", size = 7.5),
      axis.ticks = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(family = "serif", size = 9.5, hjust = 0.5, face = "bold")
    )
  if (!show_y_labels) p <- p + theme(axis.text.y = element_blank())
  p
}

p1 <- make_col("PC1", TRUE)
p2 <- make_col("PC2", FALSE) + labs(x = "canonical factor score")
p3 <- make_col("PC3", FALSE)

fig <- (p1 | p2 | p3) +
  plot_annotation(
    title = paste0("Canonical factor scores with 90% credible intervals: each target (diamond), its three ",
                    "closest style-neighbours\n(circles), and the three most distant players (grey squares)"),
    theme = theme(plot.title = element_text(size = 10, face = "bold", family = "serif"))
  )

ggsave(file.path(paths$figures, "50_score_forest_plot.png"), fig, width = 9.0, height = 8.4, dpi = 200)
message("Stage 50 complete: images/50_score_forest_plot.png")
