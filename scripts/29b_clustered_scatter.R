suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

scores <- read.csv("outputs/tables/29_player_factor_scores.csv",
                   stringsAsFactors = FALSE)

set.seed(42)
km <- kmeans(cbind(scores$pc1_score, scores$pc2_score), centers = 5, nstart = 50)
scores$cluster <- km$cluster

cents <- as.data.frame(km$centers)
colnames(cents) <- c("PC1", "PC2")
cents$cluster <- seq_len(nrow(cents))

# PC1 interpretation (confirmed from known players):
#   NEGATIVE PC1 = attacking output  (Cristiano=-7.5, Neymar=-8.3, Isak, Stuani)
#   POSITIVE PC1 = defensive work    (Ramos=+5.6, Piqué=+6.2, Busquets=+4.7)
# PC2 interpretation:
#   NEGATIVE PC2 = technical/creative (Xavi=-7.4, Iniesta=-6.8, Neymar=-6.7)
#   POSITIVE PC2 = direct/physical    (Dovbyk=+5.1, Andone=+5.5)
#
# Labels assigned by ranking clusters on PC1 (most negative = most attacking):
#   Rank 1 (lowest PC1)  = Forwards          (goal-scorers: Isak, Suárez, Stuani)
#   Rank 2               = Attacking MFs     (playmakers: Neymar, Isco; PC2 negative)
#   Rank 3               = Midfielders       (Xavi, Iniesta, Carvajal, Alba)
#   Rank 4               = Defensive MFs     (Busquets, Mascherano, Blind)
#   Rank 5 (highest PC1) = Centre backs      (Piqué, Ramos, Albiol)

pc1_rank   <- rank(cents$PC1)   # 1 = most negative = most attacking
archetypes <- c("Forwards", "Attacking midfielders", "Midfielders",
                "Defensive midfielders", "Centre backs")
label_map  <- archetypes[pc1_rank]

message("Centroid → archetype mapping:")
print(data.frame(cents[order(cents$PC1), ], archetype = archetypes))

scores$archetype <- label_map[scores$cluster]

archetype_colors <- c(
  "Forwards"              = "#E41A1C",
  "Attacking midfielders" = "#FF7F00",
  "Midfielders"           = "#4DAF4A",
  "Defensive midfielders" = "#377EB8",
  "Centre backs"          = "#984EA3"
)

# Label only the most extreme players per cluster (clear archetypes)
# Top 2 by |PC1| per cluster + top 2 by |PC2| per cluster, deduplicated, capped at 28
label_pts <- bind_rows(
  scores |> group_by(archetype) |> slice_min(order_by = pc1_score, n = 3) |> ungroup(),
  scores |> group_by(archetype) |> slice_max(order_by = pc1_score, n = 3) |> ungroup(),
  scores |> group_by(archetype) |> slice_min(order_by = pc2_score, n = 2) |> ungroup(),
  scores |> group_by(archetype) |> slice_max(order_by = pc2_score, n = 2) |> ungroup()
) |>
  filter(!is.na(player_name)) |>
  distinct(player_name, .keep_all = TRUE) |>
  mutate(extremity = pc1_score^2 + pc2_score^2) |>
  slice_max(order_by = extremity, n = 28)

message("Labelled players:")
print(label_pts[order(label_pts$pc1_score), c("player_name", "archetype", "pc1_score", "pc2_score")])

p <- ggplot(scores, aes(pc1_score, pc2_score, color = archetype)) +
  geom_point(alpha = 0.22, size = 0.85) +
  geom_text_repel(
    data = label_pts,
    aes(label = player_name),
    size = 2.5, max.overlaps = 35,
    box.padding = 0.38, point.padding = 0.2,
    min.segment.length = 0.2, seed = 42
  ) +
  scale_color_manual(values = archetype_colors, name = "Archetype") +
  labs(
    title    = "Player archetypes in factor space — stage 28 (K=3, P=48)",
    subtitle = paste0(
      "PC1 (61.8%): attacking output ← | → defensive work-rate\n",
      "PC2 (27.7%): technical / creative ← | → direct / physical"
    ),
    x = "PC1 (61.8% of common variance)",
    y = "PC2 (27.7% of common variance)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 8.5, color = "grey35",
                                    lineheight = 1.4)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5), nrow = 1))

ggsave("outputs/figures/29_factor_score_scatter_clustered.png", p,
       width = 10, height = 8, dpi = 150)
message("Saved: outputs/figures/29_factor_score_scatter_clustered.png")
