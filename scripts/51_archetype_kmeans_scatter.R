# Stage 51 -- 3-panel k-means archetype scatter (Results-chapter §7 revision).
#
# §7's prose describes a k-means clustering result (5 archetypes on the PC1-PC2 plane)
# without ever showing it -- a close version already exists in the professor deck
# (docs/professor/professor_summary.pdf's Figure 7, images/29_factor_score_scatter_
# clustered.png, from scripts/29b_clustered_scatter.R), but only as a single PC1-PC2
# panel in that deck's own casual ggplot style. This script rebuilds the same result as
# a proper 3-panel figure (PC1-PC2, PC1-PC3, PC2-PC3) in the thesis's own visual
# conventions, matching scripts/48_style_space_scatter.R's patchwork structure and theme.
#
# Clustering is refit here (same seed, same k, same variables) rather than read back from
# 29b_player_archetypes.csv, so the figure is self-contained and always matches whatever
# 29_player_factor_scores.csv currently holds -- it is fit ONCE on the PC1-PC2 plane only
# (matching what the §7 prose claims); the other two panels reuse that same per-player
# archetype label and just plot a different pair of axes, not a fresh clustering.
#
# Reads:  outputs/tables/29_player_factor_scores.csv (all 1,529 players' posterior-mean PC1-3)
# Writes: outputs/figures/51_archetype_kmeans_scatter.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/51_archetype_kmeans_scatter.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "patchwork"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork) })

scores <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"), show_col_types = FALSE) |>
  filter(!is.na(player_name))

# ---- k-means on the PC1-PC2 plane only (matches the §7 prose) ---------------------------
set.seed(42)
km <- kmeans(cbind(scores$pc1_score, scores$pc2_score), centers = 5, nstart = 50)
scores$cluster <- km$cluster

cents <- as.data.frame(km$centers)
colnames(cents) <- c("PC1", "PC2")

# Archetypes assigned by ranking cluster centroids on PC1 (most negative = most attacking),
# identical convention to scripts/29b_clustered_scatter.R.
pc1_rank   <- rank(cents$PC1)
ARCHETYPES <- c("Forwards", "Attacking midfielders", "Midfielders",
                 "Defensive midfielders", "Centre backs")
scores$archetype <- factor(ARCHETYPES[pc1_rank[scores$cluster]], levels = ARCHETYPES)

message("Centroid -> archetype mapping:")
print(data.frame(cents[order(cents$PC1), ], archetype = ARCHETYPES))

ARCH_COL <- c(
  "Forwards"              = "#E41A1C",
  "Attacking midfielders" = "#FF7F00",
  "Midfielders"           = "#4DAF4A",
  "Defensive midfielders" = "#377EB8",
  "Centre backs"          = "#984EA3"
)

# ---- three panels, shared legend, thesis theme (matches stage 48) -----------------------
make_panel <- function(a_var, b_var, a_lab, b_lab) {
  ggplot(scores, aes(x = .data[[a_var]], y = .data[[b_var]], colour = archetype)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "#9E9E9E") +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#9E9E9E") +
    geom_point(size = 0.9, alpha = 0.45, stroke = 0) +
    scale_colour_manual(values = ARCH_COL, name = "Archetype", drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.6))) +
    labs(x = a_lab, y = b_lab) +
    scale_x_continuous(expand = expansion(mult = 0.10)) +
    scale_y_continuous(expand = expansion(mult = 0.10)) +
    theme_minimal(base_size = 10) +
    theme(axis.title = element_text(family = "serif"),
          axis.text = element_text(family = "serif", size = 7.5),
          legend.text = element_text(family = "serif", size = 8.5),
          legend.title = element_text(family = "serif", size = 9, face = "bold"),
          panel.grid.minor = element_blank())
}

p1 <- make_panel("pc1_score", "pc2_score", "Defensive engagement", "Technical orientation")
p2 <- make_panel("pc1_score", "pc3_score", "Defensive engagement", "Possession retention")
p3 <- make_panel("pc2_score", "pc3_score", "Technical orientation", "Possession retention")

# Title deliberately omitted (G2): the LaTeX caption states what the panels and colours are.
fig <- (p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(file.path(paths$figures, "51_archetype_kmeans_scatter.png"), fig, width = 9.6, height = 3.7, dpi = 200)
message("Stage 51 complete: images/51_archetype_kmeans_scatter.png")
