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

# Print centroid positions so label mapping is transparent
cents <- as.data.frame(km$centers)
colnames(cents) <- c("PC1", "PC2")
cents$cluster <- seq_len(nrow(cents))
message("Cluster centroids:")
print(cents[order(cents$PC1), ])

# Assign archetype labels based on centroid position
# PC1 = offensive-defensive axis (high = attack, low = defend)
# PC2 = width / creativity (high = wide/creative)
label_map <- character(5)
for (k in 1:5) {
  p1 <- cents$PC1[k]; p2 <- cents$PC2[k]
  label_map[k] <- dplyr::case_when(
    p1 >  0.5 & p2 >  0.3 ~ "Attacking / wide",
    p1 >  0.5 & p2 <= 0.3 ~ "Forwards",
    p1 < -0.5 & p2 >  0.3 ~ "Fullbacks",
    p1 < -0.5              ~ "Defenders / CDMs",
    p2 >  0.5              ~ "Wide midfielders",
    TRUE                   ~ "Central midfielders"
  )
}
message("Archetype labels:")
print(data.frame(cluster = 1:5, label = label_map))

scores$archetype <- label_map[scores$cluster]

# Keep only representative labels: top 3 by |PC1| and top 3 by |PC2| per cluster
extreme_pc1 <- scores |>
  group_by(cluster) |>
  slice_max(order_by = abs(pc1_score), n = 3) |>
  ungroup()
extreme_pc2 <- scores |>
  group_by(cluster) |>
  slice_max(order_by = abs(pc2_score), n = 3) |>
  ungroup()

label_pts <- bind_rows(extreme_pc1, extreme_pc2) |>
  distinct(player_name, .keep_all = TRUE) |>
  # Keep only the most extreme examples overall to avoid crowding
  slice_max(order_by = pc1_score^2 + pc2_score^2, n = 30)

archetype_colors <- c(
  "Forwards"           = "#E41A1C",
  "Attacking / wide"   = "#FF7F00",
  "Wide midfielders"   = "#4DAF4A",
  "Central midfielders"= "#377EB8",
  "Defenders / CDMs"   = "#984EA3",
  "Fullbacks"          = "#A65628"
)

p <- ggplot(scores, aes(pc1_score, pc2_score, color = archetype)) +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_text_repel(
    data = label_pts,
    aes(label = player_name),
    size = 2.5, max.overlaps = 30,
    box.padding = 0.35, point.padding = 0.2,
    min.segment.length = 0.2, seed = 42
  ) +
  scale_color_manual(values = archetype_colors, drop = FALSE,
                     name = "Cluster archetype") +
  labs(
    title    = "Player archetypes in factor space (K=3 model, stage 28)",
    subtitle = "PC1: offensive volume vs. defensive work-rate  |  PC2: width / creativity",
    x = sprintf("PC1 (%.1f%% of shared variance)", 61.8),
    y = sprintf("PC2 (%.1f%% of shared variance)", 27.7)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9, color = "grey40")
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5)))

ggsave("outputs/figures/29_factor_score_scatter_clustered.png", p,
       width = 10, height = 8, dpi = 150)
message("Saved: outputs/figures/29_factor_score_scatter_clustered.png")
