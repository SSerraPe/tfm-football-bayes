# Stage 48 -- Style-space scatter (Results-chapter rewrite) + similarity table fragment
# + the "3 most distant players" selection used by stage 50's forest plot.
#
# Three panels (PC1-PC2, PC1-PC3, PC2-PC3): all 1,529 players' posterior-mean canonical
# scores as a grey population cloud, per-target style-neighbours (small coloured circles),
# and the target itself (a larger coloured diamond, the only labelled point). Replaces the
# single PC1-vs-PC2 k-means-coloured scatter (stage 29b) and the top/bottom-4-per-axis
# table's role as "the" style-space figure -- Table 8 (43_archetype_combined.tex) is kept
# unchanged alongside this, since it answers a different question (who is most extreme)
# that this figure deliberately does not (only the 3 targets are labelled here).
#
# Distant-player rule (documented, not just computed): for each target, style distance
# (Euclidean, standardised posterior-mean PC1-3, identical to stage 37's metric) to all
# other 1,528 players; among the top decile of that distance, the three most recognisable
# players, capped at not reaching past roughly the top 30 by raw distance. Applied
# identically to all three targets. For Messi specifically, no recognisable player appears
# within that range -- the entire top of his distance ranking is low-minute fringe
# defenders/wide players -- so the honest top-3 by raw distance are used instead and
# reported as such, rather than reaching further to manufacture a recognisable name. This
# mirrors the precedent already set for Table 8's PC3-top row (an earlier revision phase of
# this same thesis): an honest, slightly less glamorous answer over a cherry-picked one.
# The recognisability judgement itself is manual (there is no computable proxy for
# "well-known footballer" in this dataset -- number of seasons observed does not track it;
# checked directly, see git history for the exploration), exactly as Table 8's own
# selection was manual and disclosed.
#
# Reads:  outputs/tables/29_player_factor_scores.csv (all 1,529 players' posterior-mean
#           PC1-3, for the population cloud and the distance computation)
#         outputs/tables/37_player_similarity_top10.csv (top-10 neighbours per player,
#           reused for the 3 nearest neighbours per target)
#         outputs/tables/46_case_studies_posterior_mean.csv (level_gap, for the table)
# Writes: images/48_style_space_scatter.png
#         outputs/tables/48_distant_players.csv (target, player, rank, style_dist -- for
#           stage 49/50 and for disclosure in the thesis text)
#         tables/48_similarity_table.tex (booktabs/multirow/siunitx fragment)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/48_style_space_scatter.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "patchwork", "ggrepel"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(ggrepel) })

# SIGN CORRECTION (found during this revision): 29_outlier_study.R's eigendecomposition
# (`eig <- eigen(LLt); V <- eig$vectors[,1:K]`) is never sign-anchored, unlike stage 31's
# (`31_pca_with_ci.R`, which anchors each PC so the feature with the largest |loading| is
# positive -- the convention Remark 4.24/rem:sign specifies and the one the new loading
# heatmap, stage 47, is built on). Checked directly: stage 29's PC1 happens to already
# agree with that convention (coincidentally, whatever LAPACK returned), but its PC2 and
# PC3 are the exact reflection of it -- e.g. Xavi, an elite deep-lying passer, comes out at
# pc2_score = -7.40 under stage 29's raw sign, when a passing-heavy player should score
# strongly POSITIVE on a "passing volume/quality (+) vs. aerial/physical (-)" axis. Left
# uncorrected, this figure's population cloud (and Table 8, patched separately in
# scripts/43_results_tables.R) would visibly contradict the loading heatmap's own sign
# convention within the same section. A sign flip is an isometry, so no distance, ranking,
# or nearest-neighbour result anywhere in this thesis is affected -- only which quadrant
# things are plotted in.
scores <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"), show_col_types = FALSE) |>
  filter(!is.na(player_name)) |>
  mutate(pc2_score = -pc2_score, pc3_score = -pc3_score)
top10  <- read_csv(file.path(paths$tables, "37_player_similarity_top10.csv"), show_col_types = FALSE)
lvl    <- read_csv(file.path(paths$tables, "46_case_studies_posterior_mean.csv"), show_col_types = FALSE)

TARGETS <- c("L. Messi", "Sergio Ramos", "Xavi")
TARGET_COL <- c("L. Messi" = "#C0392B", "Sergio Ramos" = "#1F6F8B", "Xavi" = "#7B3294")
# Xavi's colour is a violet substitute for the brief's dark-goldenrod (#B8860B): under
# deuteranopia/protanopia (the most common colour-vision deficiencies), brick-red and dark
# goldenrod both lose their red/green contribution and collapse toward a similar muddy
# brown, while violet stays clearly separated from both the red and the teal-blue under
# all three common CVD simulations. Used consistently in this figure and in stage 50.

neighbours <- lapply(TARGETS, function(tg) {
  top10 |> filter(player_name == tg, rank <= 3) |>
    transmute(target = tg, player_name = similar_player, rank)
}) |> bind_rows()

# ---- distant players: computed distance ranking (for disclosure) + the manual, disclosed
#      recognisability pick described in the header comment -----------------------------
pc_scaled <- scale(as.matrix(scores[, c("pc1_score", "pc2_score", "pc3_score")]))
Dmat <- as.matrix(dist(pc_scaled, method = "euclidean"))
rownames(Dmat) <- colnames(Dmat) <- as.character(scores$player_id)
name_to_id <- setNames(as.character(scores$player_id), scores$player_name)
name_to_id <- name_to_id[!duplicated(scores$player_name)]

DISTANT_PICKS <- list(
  "L. Messi"     = c("O. Vranješ", "Juan Iglesias", "Iván Alejo"),      # raw top-3 (no recognisable player in-decile)
  "Sergio Ramos" = c("Hélder Costa", "L. Messi", "Pedro Porro"),          # recognisable, within top 15 by distance
  "Xavi"         = c("A. Budimir", "C. Stuani", "Y. En-Nesyri")           # recognisable, within top 30 by distance
)

distant_df <- lapply(TARGETS, function(tg) {
  pid <- name_to_id[[tg]]
  dvec <- Dmat[pid, ]; dvec <- dvec[names(dvec) != pid]
  picks <- DISTANT_PICKS[[tg]]
  tibble(target = tg, player_name = picks, rank = seq_along(picks),
         style_dist = round(dvec[name_to_id[picks]], 4))
}) |> bind_rows()
write_csv(distant_df, file.path(paths$tables, "48_distant_players.csv"))
message("Wrote 48_distant_players.csv:")
print(as.data.frame(distant_df))

# ---- scatter -------------------------------------------------------------------------------
pos <- scores |> select(player_name, pc1_score, pc2_score, pc3_score)

make_panel <- function(a_var, b_var, a_lab, b_lab) {
  p <- ggplot(pos, aes(x = .data[[a_var]], y = .data[[b_var]])) +
    geom_point(size = 1.0, colour = "#D2D2D2", alpha = 0.55, stroke = 0) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "#9E9E9E") +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9E9E9E")
  target_pos_all <- pos |> filter(player_name %in% TARGETS) |>
    mutate(colour = unlist(TARGET_COL[player_name]))
  for (tg in TARGETS) {
    col <- TARGET_COL[[tg]]
    nb_names <- neighbours |> filter(target == tg) |> pull(player_name)
    nb_pos <- pos |> filter(player_name %in% nb_names)
    tg_pos <- pos |> filter(player_name == tg)
    p <- p +
      geom_point(data = nb_pos, aes(x = .data[[a_var]], y = .data[[b_var]]),
                 size = 2.6, colour = col, alpha = 0.75, shape = 21, fill = col, stroke = 0.4) +
      geom_point(data = tg_pos, aes(x = .data[[a_var]], y = .data[[b_var]]),
                 size = 4.2, colour = col, shape = 23, fill = col, stroke = 0.9)
  }
  p <- p +
    geom_text_repel(data = target_pos_all, aes(x = .data[[a_var]], y = .data[[b_var]], label = player_name),
                     colour = target_pos_all$colour, fontface = "bold", size = 3.0, family = "serif",
                     min.segment.length = 0.3, segment.linewidth = 0.35, box.padding = 0.45,
                     point.padding = 0.15, seed = 1, force = 3, max.overlaps = Inf)
  p + labs(x = a_lab, y = b_lab) +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    scale_y_continuous(expand = expansion(mult = 0.10)) +
    theme_minimal(base_size = 10) +
    theme(axis.title = element_text(family = "serif"),
          axis.text = element_text(family = "serif", size = 7.5),
          panel.grid.minor = element_blank())
}

p1 <- make_panel("pc1_score", "pc2_score", "PC1", "PC2")
p2 <- make_panel("pc1_score", "pc3_score", "PC1", "PC3")
p3 <- make_panel("pc2_score", "pc3_score", "PC2", "PC3")

fig <- (p1 | p2 | p3) +
  plot_annotation(
    title = "Canonical style space: case-study players (diamonds) and their three closest style-neighbours (circles)",
    theme = theme(plot.title = element_text(size = 10.5, face = "bold", family = "serif"))
  )

ggsave(file.path(paths$figures, "48_style_space_scatter.png"), fig, width = 9.4, height = 3.4, dpi = 200)
message("Stage 48 complete: images/48_style_space_scatter.png")

# ---- similarity table fragment (booktabs + multirow + siunitx) ----------------------------
tbl_rows <- lapply(TARGETS, function(tg) {
  lvl |> filter(query_player_name == tg, rank <= 3) |> arrange(rank) |>
    transmute(target = tg, neighbor = neighbor_player_name, style_dist = style_dist, level_gap = level_gap)
}) |> bind_rows()

fmt <- function(x) sprintf("%.3f", x)

# Numeric columns use siunitx's S type for decimal alignment; target/neighbour stay plain
# text (l). \multirow marks the first row of each target's 3-row block.
tex <- c(
  "\\begin{tabular}{l l S[table-format=1.3] S[table-format=1.3]}",
  "\\toprule",
  "Target & Similar player & {Style distance} & {Level gap} \\\\",
  "\\midrule"
)
for (tg in TARGETS) {
  sub <- tbl_rows |> filter(target == tg)
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    target_cell <- if (i == 1) sprintf("\\multirow{3}{*}{%s}", tg) else ""
    tex <- c(tex, sprintf("%s & %s & %s & %s \\\\", target_cell, r$neighbor, fmt(r$style_dist), fmt(r$level_gap)))
  }
  if (tg != tail(TARGETS, 1)) tex <- c(tex, "\\midrule")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")

out_dir <- file.path(paths$scripts, "..", "tfm_mesio_latex_2026", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(tex, file.path(out_dir, "48_similarity_table.tex"))
message("Wrote tables/48_similarity_table.tex")
