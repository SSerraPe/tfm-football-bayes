# Stage 50 -- Per-target forest-plot "cards" of canonical factor scores under the fixed
# reference orientation (§8 revision, 2026-09).
#
# Previously one combined 4x3 grid on a single, globally-shared x-scale across all four
# targets (a tight cluster like Ramos's looked compressed next to Messi's much wider
# spread). Now four independent PNGs, one per target, each with its own x-scale (shared
# across that target's own PC1/PC2/PC3 panels, since canonical scores are on directly
# comparable footing -- but not shared across targets any more), embedded one per player
# inside that player's own §8.3 subsubsection rather than all four stacked in one figure.
#
# Each card: 7 rows (target + 3 style-neighbours + 3 most distant players, same
# composition and row order as before), 3 columns (PC1/PC2/PC3). The similarity score
# (stage 48) is appended directly to each non-target row's y-axis label ("Neymar (45%)"),
# shown once per row (on the leftmost, PC1 panel, where row labels already live) rather
# than repeated on all three panels.
#
# Reads:  outputs/tables/49_score_credible_intervals.csv (player x pc x {q05,q50,q95})
#         outputs/tables/48_distant_players.csv, outputs/tables/48_similarity_scores.csv
#         outputs/tables/37_player_similarity_top10.csv (row order/composition per target)
# Writes: images/50_score_forest_<slug>.png, one per target

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/50_score_forest_plot.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "patchwork"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

ci      <- read_csv(file.path(paths$tables, "49_score_credible_intervals.csv"), show_col_types = FALSE)
distant <- read_csv(file.path(paths$tables, "48_distant_players.csv"), show_col_types = FALSE)
simsc   <- read_csv(file.path(paths$tables, "48_similarity_scores.csv"), show_col_types = FALSE)
top10   <- read_csv(file.path(paths$tables, "37_player_similarity_top10.csv"), show_col_types = FALSE)

TARGETS <- c("L. Messi", "Sergio Ramos", "Xavi", "C. Stuani")
TARGET_SLUG <- c("L. Messi" = "messi", "Sergio Ramos" = "ramos", "Xavi" = "xavi", "C. Stuani" = "stuani")
TARGET_COL <- c("L. Messi" = "#C0392B", "Sergio Ramos" = "#1F6F8B", "Xavi" = "#7B3294", "C. Stuani" = "#2E7D32")
DISTANT_COL <- "#7A7A7A"

build_card <- function(tg) {
  nb <- top10 |> filter(player_name == tg, rank <= 3) |> arrange(rank) |> pull(similar_player)
  di <- distant |> filter(target == tg) |> arrange(rank) |> pull(player_name)
  row_order <- tibble(
    target = tg, player_name = c(tg, nb, di),
    role = c("target", rep("neighbour", 3), rep("distant", 3)),
    row_in_block = 1:7
  ) |>
    mutate(y = -row_in_block) |>
    left_join(simsc |> filter(target == tg) |> select(player_name, similarity_score),
               by = "player_name") |>
    mutate(row_label = ifelse(is.na(similarity_score), player_name,
                                sprintf("%s (%.0f%%)", player_name, similarity_score)))

  plot_df <- row_order |>
    left_join(ci, by = "player_name", relationship = "many-to-many") |>
    mutate(
      colour = ifelse(role == "distant", DISTANT_COL, TARGET_COL[[tg]]),
      shape  = case_when(role == "target" ~ 23, role == "neighbour" ~ 21, role == "distant" ~ 22),
      size   = case_when(role == "target" ~ 4.2, role == "neighbour" ~ 3.0, role == "distant" ~ 3.0)
    )

  # x-limits shared across this target's own 3 panels (PCs are on comparable footing),
  # but computed from this target's 7 rows only -- not shared with the other 3 targets.
  xmax <- max(abs(c(plot_df$q05, plot_df$q95))) * 1.12

  # Panel titles use the named components (Table 5/6); pcx itself stays "PC1"/"PC2"/"PC3"
  # since that is the literal key used by the upstream CSVs' `pc` column.
  PC_LABELS <- c("PC1" = "Defensive engagement", "PC2" = "Technical orientation",
                 "PC3" = "Possession retention")

  make_col <- function(pcx, show_y_labels) {
    df <- plot_df |> filter(pc == pcx)
    p <- ggplot(df, aes(y = y)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "#AAAAAA") +
      geom_segment(aes(x = q05, xend = q95, yend = y, colour = colour),
                   linewidth = 1.5, alpha = 0.85, lineend = "round") +
      geom_point(aes(x = q50, fill = colour, shape = shape, size = size),
                 colour = "white", stroke = 0.7) +
      scale_colour_identity() + scale_fill_identity() + scale_shape_identity() + scale_size_identity() +
      scale_x_continuous(limits = c(-xmax, xmax)) +
      scale_y_continuous(breaks = row_order$y,
                          labels = if (show_y_labels) row_order$row_label else NULL,
                          expand = expansion(add = 0.6)) +
      labs(x = NULL, y = NULL, title = PC_LABELS[[pcx]]) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.y = element_text(family = "serif", size = 8.2, hjust = 1),
        axis.text.x = element_text(family = "serif", size = 7.5),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(family = "serif", size = 8.6, hjust = 0.5, face = "bold")
      )
    if (!show_y_labels) p <- p + theme(axis.text.y = element_blank())
    p
  }

  p1 <- make_col("PC1", TRUE)
  p2 <- make_col("PC2", FALSE) + labs(x = "canonical factor score")
  p3 <- make_col("PC3", FALSE)

  # Title deliberately omitted (G2): the LaTeX caption states what the diamond, circles,
  # and grey squares are, and which player this card belongs to.
  p1 | p2 | p3
}

for (tg in TARGETS) {
  fig <- build_card(tg)
  out_path <- file.path(paths$figures, sprintf("50_score_forest_%s.png", TARGET_SLUG[[tg]]))
  ggsave(out_path, fig, width = 8.2, height = 2.6, dpi = 200)
  message("Wrote ", out_path)
}
message("Stage 50 complete: four per-target forest-plot cards written.")
