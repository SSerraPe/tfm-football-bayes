# Stage 37 — Player radar plots (FIFA/FC-style, percentile) + similarity metric
#
# Produces:
#   37_radar_messi_cr7.png            — Messi vs Cristiano Ronaldo
#   37_radar_messi_neymar.png         — Messi vs Neymar
#   37_radar_mbappe_cr7.png           — Mbappé vs Cristiano Ronaldo
#   37_radar_suarez_benzema.png       — Suárez vs Benzema
#   37_radar_iniesta_busquets.png     — Iniesta vs Busquets
#   37_radar_ramos_pique.png          — Ramos vs Piqué
#   37_player_scores_pct.csv          — 8-group percentile scores for all 1529 players
#   37_player_similarity_top10.csv    — top-10 nearest neighbours per player (PC space)
#   37_similarity_example.png         — Messi's top-10 on the PC scatter
#   37_similarity_examples.txt        — sanity-check text for named players

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "scripts/37_player_radar_similarity.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
})
source(file.path(model_root, "src", "loading_visualization.R"))

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load all data
# ══════════════════════════════════════════════════════════════════════════════
eff    <- read_csv(file.path(paths$tables, "14_player_effect_means.csv"),
                   show_col_types = FALSE)
scores <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"),
                   show_col_types = FALSE)
ci_tbl <- read_csv(file.path(paths$tables, "31_pca_loading_ci.csv"),
                   show_col_types = FALSE)

eff <- eff |>
  left_join(scores |> select(player_id, player_name), by = "player_id")

cat("Players with names:", sum(!is.na(eff$player_name)), "/", nrow(eff), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 2. Custom radar group mapping
#    "Shooting rates" had 0 features in the 48-feature model (both rate_* dropped).
#    Fix: split Attacking output → ATTACK (creation/box presence) + SHOOT (finishing).
# ══════════════════════════════════════════════════════════════════════════════
radar_groups <- list(
  "Attacking output"    = c("per90_touch_in_box", "per90_assists", "per90_xg_assist"),
  "Shooting rates"      = c("per90_goals", "per90_shots", "per90_xg_shot",
                             "per90_head_shots"),
  "Passing volume"      = c("per90_received_pass", "per90_passes",
                             "per90_forward_passes", "per90_back_passes",
                             "per90_lateral_passes", "per90_long_passes",
                             "per90_smart_passes", "per90_progressive_passes",
                             "per90_passes_to_final_third", "per90_key_passes",
                             "per90_through_passes", "per90_crosses",
                             "per90_shot_assists"),
  "Passing quality"     = c("rate_successful_passes",
                             "rate_successful_progressive_passes"),
  "Progressive actions" = c("per90_dribbles", "per90_progressive_run",
                             "per90_accelerations", "rate_successful_dribbles"),
  "Defending"           = c("per90_losses", "per90_own_half_losses",
                             "per90_dangerous_own_half_losses",
                             "per90_recoveries", "per90_opponent_half_recoveries",
                             "per90_dangerous_opponent_half_recoveries",
                             "per90_interceptions", "per90_clearances",
                             "per90_shots_blocked", "per90_sliding_tackles",
                             "per90_dribbles_against", "rate_defensive_duels_won"),
  "Dueling"             = c("per90_defensive_duels", "per90_offensive_duels",
                             "per90_aerial_duels", "per90_pressing_duels",
                             "per90_loose_ball_duels", "rate_offensive_duels_won",
                             "rate_aerial_duels_won"),
  "Discipline"          = c("per90_fouls", "per90_yellow_cards")
)

group_order_vec <- names(radar_groups)
spoke_labels    <- c("ATTACK", "SHOOT", "PASSING", "QUALITY",
                     "PROG",   "DEFEND", "DUEL",   "DISCIP")
names(spoke_labels) <- group_order_vec

# ══════════════════════════════════════════════════════════════════════════════
# 3. Compute loading-weighted group composites
# ══════════════════════════════════════════════════════════════════════════════
pc1_weights <- ci_tbl |>
  filter(pc == "PC1") |>
  select(feature, abs_loading = loading) |>
  mutate(abs_loading = abs(abs_loading))

feat_cols <- names(eff)[!names(eff) %in% c("player_id", "player_index", "player_name")]

# Features where higher raw value = worse outcome: negate before computing composite
# so that higher percentile always means "better" on every spoke.
#   Discipline: fewer fouls/cards = better (fouls_suffered stays positive — being fouled is good)
#   Defending:  fewer ball losses = better
negate_features <- c(
  "per90_fouls", "per90_yellow_cards",
  "per90_losses", "per90_own_half_losses", "per90_dangerous_own_half_losses"
)

compute_group_scores <- function(player_row, feat_cols, radar_groups, pc1_weights) {
  vals <- setNames(as.numeric(player_row[feat_cols]), feat_cols)
  vals[names(vals) %in% negate_features] <- -vals[names(vals) %in% negate_features]
  vapply(group_order_vec, function(g) {
    feats_g <- intersect(radar_groups[[g]], feat_cols)
    if (length(feats_g) == 0) return(NA_real_)
    v <- vals[feats_g]
    w <- pc1_weights$abs_loading[match(feats_g, pc1_weights$feature)]
    w[is.na(w)] <- 1.0
    sum(v * w, na.rm = TRUE) / sum(w[!is.na(v)], na.rm = TRUE)
  }, numeric(1))
}

cat("Computing group composites for", nrow(eff), "players...\n")
grp_mat <- t(apply(eff, 1, compute_group_scores,
                   feat_cols = feat_cols, radar_groups = radar_groups,
                   pc1_weights = pc1_weights))
colnames(grp_mat) <- group_order_vec

# Verify no empty spokes
empty_spokes <- colSums(is.na(grp_mat)) == nrow(grp_mat)
if (any(empty_spokes))
  warning("Empty spokes: ", paste(names(empty_spokes)[empty_spokes], collapse=", "))

# ══════════════════════════════════════════════════════════════════════════════
# 4. Convert to 0-100 percentile
# ══════════════════════════════════════════════════════════════════════════════
pct_mat <- apply(grp_mat, 2, function(x) {
  rank(x, ties.method = "average", na.last = "keep") / sum(!is.na(x)) * 100
})
colnames(pct_mat) <- group_order_vec

player_pct <- eff |>
  select(player_id, player_name) |>
  bind_cols(as.data.frame(pct_mat)) |>
  left_join(scores |> select(player_id, pc1_score, pc2_score, pc3_score, n_seasons),
            by = "player_id")

write_csv(player_pct, file.path(paths$tables, "37_player_scores_pct.csv"))
cat("Saved: 37_player_scores_pct.csv\n")

# Spot-check
check <- player_pct |>
  filter(player_name %in% c("L. Messi", "Neymar", "Cristiano Ronaldo")) |>
  select(player_name, all_of(group_order_vec))
cat("Spot-check percentiles:\n"); print(check)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Radar drawing helpers
# ══════════════════════════════════════════════════════════════════════════════
n_spokes  <- length(group_order_vec)
spoke_ang <- seq(0, 2 * pi, length.out = n_spokes + 1)[-(n_spokes + 1)]

make_circle <- function(r, n_pts = 120) {
  ang <- seq(0, 2 * pi, length.out = n_pts + 1)
  data.frame(x = r * sin(ang), y = r * cos(ang), r = r)
}
ref_df   <- bind_rows(lapply(c(25, 50, 75, 100), make_circle))
label_df <- data.frame(
  lx    = 118 * sin(spoke_ang),
  ly    = 118 * cos(spoke_ang),
  label = spoke_labels[group_order_vec]
)
spoke_df <- data.frame(
  x = 0, y = 0,
  xend = 100 * sin(spoke_ang), yend = 100 * cos(spoke_ang)
)

build_polygon <- function(pname, player_pct, colour) {
  row  <- filter(player_pct, player_name == pname)[1, ]
  vals <- as.numeric(row[group_order_vec])
  vals[is.na(vals)] <- 0
  ang  <- c(spoke_ang, spoke_ang[1])
  val2 <- c(vals, vals[1])
  data.frame(x = val2 * sin(ang), y = val2 * cos(ang))
}

build_val_labels <- function(pname, player_pct) {
  row  <- filter(player_pct, player_name == pname)[1, ]
  vals <- as.numeric(row[group_order_vec])
  vals[is.na(vals)] <- 0
  vr   <- vals + 9
  data.frame(
    x     = vr * sin(spoke_ang),
    y     = vr * cos(spoke_ang),
    label = as.integer(round(vals))
  )
}

make_radar_pair <- function(p1_name, p2_name, player_pct,
                             c1 = "#e8b400", c2 = "#d73027",
                             out_path = NULL) {
  poly1 <- build_polygon(p1_name, player_pct, c1)
  poly2 <- build_polygon(p2_name, player_pct, c2)
  lab1  <- build_val_labels(p1_name, player_pct)
  lab2  <- build_val_labels(p2_name, player_pct)

  p <- ggplot() +
    geom_path(data = ref_df[ref_df$r < 100, ],
              aes(x, y, group = r), colour = "#333355", linewidth = 0.4) +
    geom_path(data = ref_df[ref_df$r == 100, ],
              aes(x, y), colour = "#4444aa", linewidth = 0.7) +
    geom_segment(data = spoke_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "#333355", linewidth = 0.5) +
    geom_polygon(data = poly2, aes(x, y),
                 fill = c2, alpha = 0.35, colour = c2, linewidth = 1.0) +
    geom_polygon(data = poly1, aes(x, y),
                 fill = c1, alpha = 0.35, colour = c1, linewidth = 1.0) +
    geom_text(data = lab1, aes(x, y, label = label),
              colour = c1, size = 3.0, fontface = "bold") +
    geom_text(data = lab2, aes(x, y, label = label),
              colour = c2, size = 3.0, fontface = "bold") +
    geom_text(data = label_df, aes(lx, ly, label = label),
              colour = "white", size = 3.2, fontface = "bold") +
    annotate("text", x = 0, y = 52, label = "50",
             colour = "#7777aa", size = 2.5) +
    # Legend via off-screen dummy points
    geom_point(data = data.frame(x = 1e6, y = 1e6,
                                 player = c(p1_name, p2_name)),
               aes(x = x, y = y, colour = player)) +
    scale_colour_manual(values = setNames(c(c1, c2), c(p1_name, p2_name))) +
    coord_fixed(xlim = c(-140, 140), ylim = c(-140, 140)) +
    labs(
      title    = sprintf("%s  vs  %s", p1_name, p2_name),
      subtitle = "Percentile among 1,529 La Liga outfield players (0 = lowest, 100 = highest)"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "#1a1a2e", colour = NA),
      panel.background = element_rect(fill = "#1a1a2e", colour = NA),
      plot.title    = element_text(colour = "white", face = "bold",
                                   hjust = 0.5, size = 14),
      plot.subtitle = element_text(colour = "#aaaacc", hjust = 0.5, size = 10),
      legend.text   = element_text(colour = "white", size = 10),
      legend.title  = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(12, 12, 12, 12)
    ) +
    guides(colour = guide_legend(
      override.aes = list(shape = 15, size = 4),
      keywidth = unit(1.2, "cm")
    ))

  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 7, height = 7.5, dpi = 150, bg = "#1a1a2e")
    cat("Saved:", out_path, "\n")
  }
  invisible(p)
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Produce all comparison pairs
# ══════════════════════════════════════════════════════════════════════════════
cat("\n--- Producing radar comparisons ---\n")

make_radar_pair("L. Messi", "Cristiano Ronaldo", player_pct,
  c1 = "#e8b400", c2 = "#c8102e",
  out_path = file.path(paths$figures, "37_radar_messi_cr7.png"))

make_radar_pair("L. Messi", "Neymar", player_pct,
  c1 = "#e8b400", c2 = "#009c3b",
  out_path = file.path(paths$figures, "37_radar_messi_neymar.png"))

make_radar_pair("K. Mbappé", "Cristiano Ronaldo", player_pct,
  c1 = "#003189", c2 = "#c8102e",
  out_path = file.path(paths$figures, "37_radar_mbappe_cr7.png"))

make_radar_pair("L. Suárez", "K. Benzema", player_pct,
  c1 = "#75aadb", c2 = "#ffffff",
  out_path = file.path(paths$figures, "37_radar_suarez_benzema.png"))

make_radar_pair("Andrés Iniesta", "Sergio Busquets", player_pct,
  c1 = "#0057a8", c2 = "#a50044",
  out_path = file.path(paths$figures, "37_radar_iniesta_busquets.png"))

make_radar_pair("Sergio Ramos", "Gerard Piqué", player_pct,
  c1 = "#ffffff", c2 = "#004d98",
  out_path = file.path(paths$figures, "37_radar_ramos_pique.png"))

# ══════════════════════════════════════════════════════════════════════════════
# 7. Similarity metric — Euclidean distance in standardised PC space
# ══════════════════════════════════════════════════════════════════════════════
cat("\n--- Computing similarity ---\n")

pc_mat <- scores |>
  filter(!is.na(player_name)) |>
  select(player_id, player_name, pc1_score, pc2_score, pc3_score) |>
  as.data.frame()

pc_scaled <- scale(as.matrix(pc_mat[, c("pc1_score", "pc2_score", "pc3_score")]))

D <- as.matrix(dist(pc_scaled, method = "euclidean"))
rownames(D) <- colnames(D) <- as.character(pc_mat$player_id)

id_to_name  <- setNames(pc_mat$player_name, as.character(pc_mat$player_id))
name_to_id  <- setNames(as.character(pc_mat$player_id), pc_mat$player_name)
name_to_id  <- name_to_id[!duplicated(pc_mat$player_name)]

find_similar_by_id <- function(pid_str, k = 10) {
  if (!pid_str %in% rownames(D)) return(NULL)
  d <- sort(D[pid_str, ])[-1]
  data.frame(
    player_id      = pid_str,
    player_name    = id_to_name[pid_str],
    similar_id     = names(head(d, k)),
    similar_player = id_to_name[names(head(d, k))],
    rank           = seq_len(k),
    distance       = round(head(d, k), 4),
    stringsAsFactors = FALSE
  )
}

find_similar <- function(player_name_str, k = 10) {
  pid <- name_to_id[player_name_str]
  if (is.na(pid) || !pid %in% rownames(D)) {
    warning("Player not found: ", player_name_str); return(NULL)
  }
  find_similar_by_id(pid, k)
}

cat("Building full similarity table...\n")
sim_df <- bind_rows(lapply(rownames(D), find_similar_by_id, k = 10))

write_csv(sim_df, file.path(paths$tables, "37_player_similarity_top10.csv"))
cat("Saved: 37_player_similarity_top10.csv\n")

# ══════════════════════════════════════════════════════════════════════════════
# 8. Similarity text examples
# ══════════════════════════════════════════════════════════════════════════════
example_players <- c("L. Messi", "Sergio Busquets", "Sergio Ramos", "K. Benzema")

lines <- c(
  "=== Stage 37: Player similarity — top 10 nearest neighbours (PC1/PC2/PC3 space) ===",
  sprintf("Date: %s", Sys.Date()), "",
  "Distance = Euclidean in standardised (PC1, PC2, PC3) space.", ""
)
for (pname in example_players) {
  sim_p <- find_similar(pname, k = 10)
  if (is.null(sim_p)) next
  lines <- c(lines,
             sprintf("--- Top 10 similar to %s ---", pname),
             sprintf("  %2d. %-35s  dist = %.3f",
                     sim_p$rank, sim_p$similar_player, sim_p$distance), "")
}
out_note <- file.path(paths$notes, "37_similarity_examples.txt")
writeLines(lines, out_note)
cat("Saved:", out_note, "\n")
cat(paste(lines, collapse = "\n"), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# 9. Similarity visualisation — Messi's top-10 on the PC scatter
# ══════════════════════════════════════════════════════════════════════════════
target <- "L. Messi"
top10  <- find_similar(target, k = 10)

scatter_df <- scores |>
  mutate(group_type = case_when(
    player_name == target                          ~ "target",
    player_name %in% top10$similar_player         ~ "top-10 similar",
    TRUE                                           ~ "other"
  ))

label_df_sim <- scatter_df |>
  filter(group_type != "other") |>
  mutate(short_name = ifelse(nchar(player_name) > 15,
                             paste0(substr(player_name, 1, 14), "."),
                             player_name))

label_colours <- ifelse(label_df_sim$group_type == "target", "#e8b400", "#d73027")

p_sim <- ggplot(scatter_df, aes(x = pc1_score, y = pc2_score)) +
  geom_point(data = ~ filter(., group_type == "other"),
             colour = "grey70", alpha = 0.25, size = 0.8) +
  geom_point(data = ~ filter(., group_type == "top-10 similar"),
             colour = "#d73027", size = 3, alpha = 0.9) +
  geom_point(data = ~ filter(., group_type == "target"),
             colour = "#e8b400", size = 5, shape = 18) +
  ggrepel::geom_text_repel(
    data = label_df_sim, aes(label = short_name),
    size = 3.0, max.overlaps = 20, seed = 42,
    colour = label_colours, box.padding = 0.35
  ) +
  labs(
    title    = sprintf("Top 10 most similar players to %s", target),
    subtitle = "Red = nearest neighbours in (PC1, PC2, PC3) archetype space. Gold = target.",
    x = "PC1 (Defensive work-rate →)",
    y = "PC2 (Technical/pass-based →)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(paths$figures, "37_similarity_example.png"),
       p_sim, width = 9, height = 6.5, dpi = 150)
cat("Saved: 37_similarity_example.png\nDone.\n")
