# SUPERSEDED (2026-09-06): this is a cleaned-up, standalone port of the original
# exploratory analysis in tfm_playground/similarity_measures.Rmd, which built the
# n=4256, p=23 player-season feature matrix and dissimilarity matrices originally
# cited in Chapters 2-3 of the thesis (tfm_latex/data_source.tex, background.tex).
#
# It has been superseded by scripts/02b_descriptive_pca_fa_cmds.R, which reruns the
# same descriptive methodology (PCA / FA / CMDS) on the corrected, final production
# dataset (N=4586, P=48, model_objects.rds) and is now the source for Chapters 2-3.
# This script is preserved for provenance only: it reproduces the earlier n=4256/p=23
# artifact (tfm_playground/datadistance_obj.Rds) from the raw data, standalone.
#
# Two bugs present in the original notebook are fixed here:
#   1. The original wrote to OUTPUTS/distance_obj.Rds (a "data/" subdirectory that was
#      never created) but the artifact actually consumed by playground_BMDS.Rmd lived
#      at tfm_playground/datadistance_obj.Rds (a different name, a different directory,
#      only reachable because it had been saved by hand outside the visible chunk
#      history). This script writes to a single, consistent path.
#   2. The original notebook's own comment ("traemos la data desde temporada 14/15
#      hasta 24/25") does not match what its code does: no season_id filter is ever
#      applied, so all 12 raw seasons (2011/12-2022/23) are used. This script keeps
#      that behaviour (all seasons) and documents it accurately instead of repeating
#      the stale claim.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/archive/exploratory_p23_similarity_analysis.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c("dplyr", "stringr"))

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

df <- read.csv(paths$raw)

# ---------------------------------------------------------------------------
# 1) Goalkeeper exclusion (match-level activity markers, as in the SQL query)
# ---------------------------------------------------------------------------
df2 <- df |>
  mutate(
    is_gk = (coalesce(total_gk_shots_against, 0) > 0) |
      (coalesce(total_gk_saves, 0) > 0) |
      (coalesce(total_gk_conceded_goals, 0) > 0) |
      (coalesce(total_gk_clean_sheets, 0) > 0) |
      (coalesce(total_gk_exits, 0) > 0)
  ) |>
  filter(!is_gk)

# ---------------------------------------------------------------------------
# 2) Aggregate match-player rows to player-season totals; 600-minute exposure filter
# ---------------------------------------------------------------------------
ps <- df2 |>
  group_by(player_id, player_name, season_id, competition_id) |>
  summarise(
    minutes = sum(total_minutes_on_field, na.rm = TRUE),
    matches = n_distinct(match_id),
    starts  = sum(coalesce(total_matches_in_start, 0), na.rm = TRUE),
    goals   = sum(coalesce(total_goals, 0), na.rm = TRUE),
    assists = sum(coalesce(total_assists, 0), na.rm = TRUE),
    shots   = sum(coalesce(total_shots, 0), na.rm = TRUE),
    xg_shot = sum(coalesce(total_xg_shot, 0), na.rm = TRUE),
    xg_ast  = sum(coalesce(total_xg_assist, 0), na.rm = TRUE),
    passes = sum(coalesce(total_passes, 0), na.rm = TRUE),
    suc_passes = sum(coalesce(total_successful_passes, 0), na.rm = TRUE),
    prog_passes = sum(coalesce(total_progressive_passes, 0), na.rm = TRUE),
    suc_prog_passes = sum(coalesce(total_successful_progressive_passes, 0), na.rm = TRUE),
    dribbles = sum(coalesce(total_dribbles, 0), na.rm = TRUE),
    suc_dribbles = sum(coalesce(total_successful_dribbles, 0), na.rm = TRUE),
    def_actions = sum(coalesce(total_defensive_actions, 0), na.rm = TRUE),
    suc_def_actions = sum(coalesce(total_successful_defensive_actions, 0), na.rm = TRUE),
    duels = sum(coalesce(total_duels, 0), na.rm = TRUE),
    duels_won = sum(coalesce(total_duels_won, 0), na.rm = TRUE),
    aerial = sum(coalesce(total_aerial_duels, 0), na.rm = TRUE),
    aerial_won = sum(coalesce(total_aerial_duels_won, 0), na.rm = TRUE),
    press_duels = sum(coalesce(total_pressing_duels, 0), na.rm = TRUE),
    press_won   = sum(coalesce(total_pressing_duels_won, 0), na.rm = TRUE),
    recov = sum(coalesce(total_recoveries, 0), na.rm = TRUE),
    opp_half_recov = sum(coalesce(total_opponent_half_recoveries, 0), na.rm = TRUE),
    touches_box = sum(coalesce(total_touch_in_box, 0), na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(minutes >= 600)

# ---------------------------------------------------------------------------
# 3) Feature engineering: per-90 volume + efficiency ratios
# ---------------------------------------------------------------------------
safe_div <- function(a, b) ifelse(b > 0, a / b, NA_real_)
per90 <- function(x, minutes) 90 * x / minutes

ps_feat <- ps |>
  mutate(
    goals_p90 = per90(goals, minutes), assists_p90 = per90(assists, minutes),
    shots_p90 = per90(shots, minutes), xg_p90 = per90(xg_shot, minutes),
    xgA_p90 = per90(xg_ast, minutes), passes_p90 = per90(passes, minutes),
    prog_passes_p90 = per90(prog_passes, minutes), dribbles_p90 = per90(dribbles, minutes),
    def_actions_p90 = per90(def_actions, minutes), duels_p90 = per90(duels, minutes),
    aerial_p90 = per90(aerial, minutes), press_duels_p90 = per90(press_duels, minutes),
    recov_p90 = per90(recov, minutes), opp_half_recov_p90 = per90(opp_half_recov, minutes),
    touches_box_p90 = per90(touches_box, minutes),
    pass_acc = safe_div(suc_passes, passes), prog_pass_acc = safe_div(suc_prog_passes, prog_passes),
    dribble_succ = safe_div(suc_dribbles, dribbles), def_succ = safe_div(suc_def_actions, def_actions),
    duel_win = safe_div(duels_won, duels), aerial_win = safe_div(aerial_won, aerial),
    press_win = safe_div(press_won, press_duels), xg_per_shot = safe_div(xg_shot, shots)
  )

# ---------------------------------------------------------------------------
# 4) Football blocks, log1p on volume features, Z-scoring, equal block weighting
# ---------------------------------------------------------------------------
blocks <- list(
  offense = c("xg_p90", "goals_p90", "assists_p90", "shots_p90", "xgA_p90", "touches_box_p90", "xg_per_shot"),
  progression_creation = c("prog_passes_p90", "passes_p90", "prog_pass_acc", "pass_acc"),
  carrying_1v1 = c("dribbles_p90", "dribble_succ"),
  defense_press = c("def_actions_p90", "def_succ", "recov_p90", "opp_half_recov_p90", "press_duels_p90", "press_win"),
  duels_phys = c("duels_p90", "duel_win", "aerial_p90", "aerial_win")
)
all_feats <- unique(unlist(blocks))

scale_z <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

id_cols <- c("player_id", "player_name", "season_id", "competition_id", "minutes", "matches")
log1p_feats <- c(
  "xg_p90", "goals_p90", "assists_p90", "shots_p90", "xgA_p90", "touches_box_p90",
  "prog_passes_p90", "passes_p90", "dribbles_p90", "def_actions_p90", "recov_p90",
  "opp_half_recov_p90", "press_duels_p90", "duels_p90", "aerial_p90"
)

X <- ps_feat |>
  select(any_of(id_cols), any_of(all_feats)) |>
  mutate(across(any_of(log1p_feats), ~ log1p(.x))) |>
  mutate(across(any_of(all_feats), scale_z))

X_mat <- X |>
  select(any_of(all_feats)) |>
  mutate(across(everything(), ~ ifelse(is.na(.x), 0, .x))) |>
  as.matrix()

block_weights <- setNames(rep(1.0, length(blocks)), names(blocks))
col_w <- setNames(rep(1, length(all_feats)), all_feats)
for (b in names(blocks)) {
  feats_b <- intersect(blocks[[b]], all_feats)
  col_w[feats_b] <- block_weights[[b]] / sqrt(length(feats_b))
}
X_w <- sweep(X_mat, 2, col_w, `*`)

# ---------------------------------------------------------------------------
# 5) Distance matrices + save artifact (single consistent path; fixes bug 1 above)
# ---------------------------------------------------------------------------
cosine_dist <- function(M) {
  M <- as.matrix(M); M[is.na(M)] <- 0
  norms <- sqrt(rowSums(M^2)); norms[norms == 0] <- 1
  Mnorm <- M / norms
  D <- 1 - (Mnorm %*% t(Mnorm))
  diag(D) <- 0
  as.dist(D)
}

dist_obj <- list(
  D_euclid = dist(X_w, method = "euclidean"),
  D_cosine = cosine_dist(X_w),
  row_keys = X |> select(any_of(id_cols)),
  X = X, X_mat = X_mat, X_w = X_w, blocks = blocks, col_w = col_w
)

out_path <- file.path(paths$processed, "archive_p23_distance_obj.rds")
saveRDS(dist_obj, out_path)
message("Wrote ", out_path, ": N=", nrow(dist_obj$X), ", P=", ncol(dist_obj$X_w),
        ", seasons=", n_distinct(dist_obj$row_keys$season_id))
