# Stage 46 — Posterior-mean-only style/level scatter for the Results chapter
# (thesis revision Phase 4: player similarity rewritten around posterior means
# only, replacing the per-draw treatment of stages 38/42/44).
#
# Style distance: Euclidean distance between two players' posterior-mean
# canonical factor scores (Definition~pca_postproc), i.e. the same K=3
# coordinates already used for the archetype clustering (stage 29b) and for
# stage 37's "classical point-estimate" nearest-neighbour table -- no per-draw
# standardisation, no credible interval on the distance itself.
#
# Level gap: Euclidean distance between two players' posterior-mean per90
# goals and per90 xG player effects (a single number per pair, not a per-draw
# posterior with a stability statistic).
#
# Reads:  outputs/tables/37_player_similarity_top10.csv (posterior-mean style
#           distance, already computed at stage 37 -- no new PC-space work)
#         outputs/tables/29_player_effect_means_k3.csv (posterior-mean player
#           effects, for per90_goals / per90_xg_shot)
#         outputs/tables/29_player_factor_scores.csv (player_index<->player_id)
#         outputs/tables/29b_player_archetypes.csv (archetype cluster labels)
# Writes: outputs/tables/46_case_studies_posterior_mean.csv
#         outputs/figures/46_style_vs_level_scatter_pm.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/46_similarity_posterior_mean.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "ggrepel"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(ggplot2) })

top10   <- read_csv(file.path(paths$tables, "37_player_similarity_top10.csv"), show_col_types = FALSE)
eff     <- read_csv(file.path(paths$tables, "29_player_effect_means_k3.csv"), show_col_types = FALSE)
idmap   <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"), show_col_types = FALSE) |>
  select(player_index, player_id)
arch    <- read_csv(file.path(paths$tables, "29b_player_archetypes.csv"), show_col_types = FALSE) |>
  select(player_id, archetype)

level_tbl <- eff |>
  left_join(idmap, by = "player_index") |>
  select(player_id, per90_goals, per90_xg_shot)

# "Xavi" added for the Results-chapter rewrite (three-target similarity treatment:
# Messi/Ramos/Xavi) -- previously absent, so the similarity table's level-gap column
# had no data for him.
CASE_STUDY_PLAYERS <- c(
  "L. Messi", "Cristiano Ronaldo", "Neymar", "L. Suárez", "K. Benzema",
  "A. Griezmann", "Iago Aspas", "Sergio Busquets", "Saúl Ñíguez",
  "L. Modrić", "D. Carvajal", "Koke", "Sergio Ramos", "Gerard Piqué",
  "J. Giménez", "Í. Martínez", "Iñaki Williams", "Xavi"
)

case_pairs <- top10 |>
  filter(player_name %in% CASE_STUDY_PLAYERS, rank <= 10) |>
  left_join(arch, by = "player_id") |>
  rename(query_player_name = player_name, query_archetype = archetype,
         neighbor_player_name = similar_player, style_dist = distance) |>
  left_join(level_tbl |> rename(player_id_q = player_id, goals_q = per90_goals, xg_q = per90_xg_shot),
             by = c("player_id" = "player_id_q")) |>
  left_join(level_tbl |> rename(player_id_n = player_id, goals_n = per90_goals, xg_n = per90_xg_shot),
             by = c("similar_id" = "player_id_n")) |>
  mutate(level_gap = sqrt((goals_q - goals_n)^2 + (xg_q - xg_n)^2)) |>
  select(query_player_name, query_archetype, neighbor_player_name, rank, style_dist, level_gap)

write_csv(case_pairs, file.path(paths$tables, "46_case_studies_posterior_mean.csv"))
message("Saved: 46_case_studies_posterior_mean.csv (", nrow(case_pairs), " rows)")

# Quick console check for the two named case studies discussed in prose
message("\n--- Benzema top neighbours (posterior-mean) ---")
print(case_pairs |> filter(query_player_name == "K. Benzema") |> arrange(rank), n = 5)
message("\n--- Ramos top neighbours (posterior-mean) ---")
print(case_pairs |> filter(query_player_name == "Sergio Ramos") |> arrange(rank), n = 5)

p <- ggplot(case_pairs, aes(x = style_dist, y = level_gap, color = query_archetype)) +
  geom_point(size = 1.8, alpha = 0.85) +
  labs(
    x = "Style distance (posterior-mean canonical factor scores)",
    y = "Level gap: |goals, xG per90| (posterior mean)",
    color = "Query archetype",
    title = "Style distance vs. level gap, top-10 style-neighbours of 17 case-study players"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(paths$figures, "46_style_vs_level_scatter_pm.png"), p, width = 8.5, height = 6, dpi = 150)
message("Stage 46 complete: outputs/figures/46_style_vs_level_scatter_pm.png")
