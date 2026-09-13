# Stage 42 — Extended Bayesian player-similarity analysis (deep dive for Results §5.8).
#
# Builds on stage 38's cached per-draw style scores (fits/38_style_score_draws.rds --
# reused as-is, avoiding another ~20-minute CmdStan-CSV column-pass) to:
#
#   (a) Broaden the named case-study set from 7 to ~17 players, spanning all five
#       archetype clusters (per 29b_player_archetypes.csv), not just attackers.
#   (b) A population-level reliability check: for a random sample of players (not
#       just named stars), what fraction have a *stable* top-1 style-neighbour
#       across posterior draws? Uses a mean+argmin-count accumulator (not full
#       per-draw storage) so it stays cheap even for ~150 query players.
#
# Kept separate from the already-committed scripts/38_bayesian_player_similarity.R
# to avoid touching validated, already-written-up work.
#
# Reads:  fits/38_style_score_draws.rds        (draws x I x K style scores, cached)
#         outputs/tables/29_player_effect_means_k3.csv  (K=3 A[i,goal/xg] source, via
#                                                          a fresh targeted column-pass)
#         outputs/tables/29b_player_archetypes.csv       (archetype cluster labels)
# Writes: outputs/tables/42_similarity_case_studies_extended.csv
#         outputs/tables/42_style_neighbor_stability_extended.csv
#         outputs/tables/42_population_reliability.csv
#         outputs/figures/42_style_vs_level_scatter_extended.png
#         outputs/figures/42_population_reliability_hist.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/42_bayesian_similarity_extended.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(ggplot2)
})

MODEL_ID <- "28_real_lowrank_a_diag_b_t_k3"
mo            <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names
P             <- length(feature_names)
GOAL_IDX      <- match("per90_goals",   feature_names)
XG_IDX        <- match("per90_xg_shot", feature_names)

F_draws <- readRDS(file.path(paths$fits, "38_style_score_draws.rds"))  # n_draws x I x K
n_draws <- dim(F_draws)[1]; I_players <- dim(F_draws)[2]; K <- dim(F_draws)[3]
message(sprintf("Loaded cached style scores: %d draws x %d players x %d factors", n_draws, I_players, K))

archetypes <- read_csv(file.path(paths$tables, "29b_player_archetypes.csv"), show_col_types = FALSE)

# ── Named case-study players: original 7 (already validated) + 10 new, spanning ──
# ── all five archetype clusters (indices confirmed against 29b_player_archetypes.csv) ──
case_study_players <- tribble(
  ~player_index, ~player_name,        ~archetype_hint,
  63,  "L. Messi",             "Attacking midfielders",
  39,  "Cristiano Ronaldo",    "Forwards",
  583, "Neymar",               "Attacking midfielders",
  33,  "Sergio Ramos",         "Centre backs",
  51,  "Gerard Piqué",         "Centre backs",
  431, "L. Suárez",            "Forwards",
  38,  "K. Benzema",           "Attacking midfielders",
  53,  "Sergio Busquets",      "Defensive midfielders",
  313, "Saúl Ñíguez",          "Defensive midfielders",
  436, "L. Modrić",            "Midfielders",
  316, "Dani Carvajal",        "Midfielders",
  99,  "Koke",                 "Midfielders",
  182, "A. Griezmann",         "Attacking midfielders",
  220, "Iago Aspas",           "Attacking midfielders",
  834, "Iñaki Williams",       "Forwards",
  807, "J. Giménez",           "Centre backs",
  173, "Íñigo Martínez",       "Centre backs"
)
stopifnot(all(case_study_players$player_index %in% archetypes$player_index))
message(sprintf("Case-study set: %d players across %d archetype clusters.",
                 nrow(case_study_players), length(unique(case_study_players$archetype_hint))))

# ── (a) Style search for the extended case-study set (full per-draw storage, ─────
# ──     one query at a time, to bound memory -- same method as stage 38) ─────────

style_neighbors <- list()
for (r in seq_len(nrow(case_study_players))) {
  qi <- case_study_players$player_index[r]; qn <- case_study_players$player_name[r]
  message("Style search for ", qn, " (player_index ", qi, ") ...")

  d_style <- matrix(NA_real_, n_draws, I_players)
  for (s in seq_len(n_draws)) {
    Fs <- F_draws[s, , ]
    Fs_std <- scale(Fs)
    d_style[s, ] <- sqrt(rowSums(sweep(Fs_std, 2, Fs_std[qi, ], `-`)^2))
  }
  d_style[, qi] <- NA_real_

  mean_d <- colMeans(d_style, na.rm = TRUE)
  q05_d  <- apply(d_style, 2, quantile, probs = 0.05, na.rm = TRUE)
  q95_d  <- apply(d_style, 2, quantile, probs = 0.95, na.rm = TRUE)

  top10_idx <- order(mean_d)[seq_len(10)]
  rank1_counts <- table(factor(apply(d_style, 1, which.min), levels = seq_len(I_players)))
  stability <- as.numeric(rank1_counts[top10_idx]) / n_draws

  style_neighbors[[qn]] <- tibble(
    query_player_index = qi, query_player_name = qn,
    query_archetype = case_study_players$archetype_hint[r],
    neighbor_index = top10_idx,
    style_dist_mean = mean_d[top10_idx], style_dist_q05 = q05_d[top10_idx], style_dist_q95 = q95_d[top10_idx],
    p_closest_neighbor = stability
  )
}
style_neighbors_tbl <- bind_rows(style_neighbors) |>
  left_join(archetypes |> select(neighbor_index = player_index, player_name, archetype),
            by = "neighbor_index")

write_csv(style_neighbors_tbl, file.path(paths$tables, "42_style_neighbor_stability_extended.csv"))
message("Saved: 42_style_neighbor_stability_extended.csv")

# ── Pass-2 A[i, goal/xg] extraction for this extended universe ───────────────

level_players <- sort(unique(c(case_study_players$player_index, style_neighbors_tbl$neighbor_index)))
message("Pass 2: extracting A[i, goals/xg_shot] for ", length(level_players), " players ...")

csv_files <- discover_cmdstan_csv_files(MODEL_ID)
py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")
find_col_indices <- function(csv_path, param_names) {
  con <- file(csv_path, "r"); on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0L) return(setNames(rep(NA_integer_, length(param_names)), param_names))
    if (!startsWith(line, "#")) {
      hdr <- strsplit(line, ",")[[1]]
      return(setNames(match(param_names, hdr), param_names))
    }
  }
}
extract_batch <- function(csv_path, chain_id, params_0idx, out_path) {
  args_pairs <- paste0(names(params_0idx), ":", params_0idx - 1L)
  cmd_args <- c(shQuote(py_script), shQuote(csv_path), as.character(chain_id),
                shQuote(out_path), args_pairs)
  ret <- system2("python3", args = cmd_args, stdout = FALSE, stderr = FALSE)
  if (ret != 0) stop("Python extraction failed for chain ", chain_id)
  read.csv(out_path, check.names = FALSE)
}
SCRATCHPAD <- file.path(tempdir(), "stage42_extract")
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

a_names <- c(paste0("A.", level_players, ".", GOAL_IDX), paste0("A.", level_players, ".", XG_IDX))
col_idx2 <- find_col_indices(csv_files[1], a_names)
missing2 <- names(col_idx2)[is.na(col_idx2)]
if (length(missing2) > 0) stop("Not found in header: ", paste(head(missing2, 5), collapse = ", "))

pass2_dfs <- lapply(seq_along(csv_files), function(i) {
  out <- file.path(SCRATCHPAD, sprintf("pass2_chain%d.csv", i))
  message("  chain ", i, " ...")
  extract_batch(csv_files[i], i, col_idx2, out)
})
pass2_df <- do.call(rbind, pass2_dfs)
goal_col <- function(i) paste0("A.", i, ".", GOAL_IDX)
xg_col   <- function(i) paste0("A.", i, ".", XG_IDX)

case_rows <- list()
for (r in seq_len(nrow(case_study_players))) {
  qi <- case_study_players$player_index[r]; qn <- case_study_players$player_name[r]
  nbrs <- style_neighbors_tbl |> filter(query_player_index == qi)
  q_goal <- pass2_df[[goal_col(qi)]]; q_xg <- pass2_df[[xg_col(qi)]]
  for (k in seq_len(nrow(nbrs))) {
    ni <- nbrs$neighbor_index[k]
    n_goal <- pass2_df[[goal_col(ni)]]; n_xg <- pass2_df[[xg_col(ni)]]
    level_gap <- sqrt((q_goal - n_goal)^2 + (q_xg - n_xg)^2)
    case_rows[[length(case_rows) + 1]] <- tibble(
      query_player_name = qn, query_archetype = case_study_players$archetype_hint[r],
      neighbor_player_name = nbrs$player_name[k], neighbor_archetype = nbrs$archetype[k],
      style_dist_mean = nbrs$style_dist_mean[k], style_dist_q05 = nbrs$style_dist_q05[k], style_dist_q95 = nbrs$style_dist_q95[k],
      p_closest_neighbor = nbrs$p_closest_neighbor[k],
      level_gap_mean = mean(level_gap), level_gap_q05 = quantile(level_gap, 0.05), level_gap_q95 = quantile(level_gap, 0.95)
    )
  }
}
case_studies_tbl <- bind_rows(case_rows) |> arrange(query_player_name, style_dist_mean)
write_csv(case_studies_tbl, file.path(paths$tables, "42_similarity_case_studies_extended.csv"))
message("Saved: 42_similarity_case_studies_extended.csv")

p <- ggplot(case_studies_tbl, aes(x = style_dist_mean, y = level_gap_mean, color = query_archetype)) +
  geom_point(size = 1.8, alpha = 0.85) +
  labs(
    title = "Style vs. level: top-10 style-neighbours of 17 case-study players across all 5 archetypes",
    subtitle = "Bottom-left = near-identical style AND output. Bottom-right = same style, very different output level.",
    x = "Style distance (posterior mean)", y = "Level gap: |goals, xG per90| (posterior mean)",
    color = "Query archetype"
  ) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(paths$figures, "42_style_vs_level_scatter_extended.png"), p, width = 9, height = 7, dpi = 150)
message("Saved: 42_style_vs_level_scatter_extended.png")

# ── (b) Population-level reliability: ~150 random players, mean+argmin-count only ──

set.seed(20260913L)
N_SAMPLE <- 150L
sample_idx <- sample(seq_len(I_players), N_SAMPLE)

dist_sum      <- matrix(0, N_SAMPLE, I_players)
closest_count <- matrix(0L, N_SAMPLE, I_players)

message("Population reliability pass: ", N_SAMPLE, " random query players x ", n_draws, " draws ...")
for (s in seq_len(n_draws)) {
  Fs <- F_draws[s, , ]
  Fs_std <- scale(Fs)
  for (qi_pos in seq_len(N_SAMPLE)) {
    qi <- sample_idx[qi_pos]
    d <- sqrt(rowSums(sweep(Fs_std, 2, Fs_std[qi, ], `-`)^2))
    d[qi] <- Inf
    dist_sum[qi_pos, ] <- dist_sum[qi_pos, ] + d
    j_star <- which.min(d)
    closest_count[qi_pos, j_star] <- closest_count[qi_pos, j_star] + 1L
  }
  if (s %% 1000 == 0) message("  draw ", s, "/", n_draws)
}
dist_mean <- dist_sum / n_draws

player_names <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"), show_col_types = FALSE) |>
  select(player_index, player_name)

reliability_rows <- lapply(seq_len(N_SAMPLE), function(qi_pos) {
  qi <- sample_idx[qi_pos]
  j_star <- which.min(dist_mean[qi_pos, ])
  tibble(
    player_index = qi,
    player_name = player_names$player_name[player_names$player_index == qi],
    top1_neighbor_index = j_star,
    top1_neighbor_name = player_names$player_name[player_names$player_index == j_star],
    top1_mean_dist = dist_mean[qi_pos, j_star],
    p_closest_neighbor = closest_count[qi_pos, j_star] / n_draws
  )
})
reliability_tbl <- bind_rows(reliability_rows)
write_csv(reliability_tbl, file.path(paths$tables, "42_population_reliability.csv"))
message("Saved: 42_population_reliability.csv")
message(sprintf("Population reliability summary: median P(closest)=%.2f, %% with P(closest)>=0.5: %.1f%%",
                 median(reliability_tbl$p_closest_neighbor),
                 100 * mean(reliability_tbl$p_closest_neighbor >= 0.5)))

p_hist <- ggplot(reliability_tbl, aes(p_closest_neighbor)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.8) +
  labs(
    title = sprintf("Stability of the top-1 style-neighbour across %d posterior draws", n_draws),
    subtitle = sprintf("Random sample of %d players (not just named stars). Median P(closest) = %.2f.",
                        N_SAMPLE, median(reliability_tbl$p_closest_neighbor)),
    x = "P(this player's mean-based top-1 neighbour is the single closest, per draw)", y = "Count"
  ) +
  theme_minimal(base_size = 11)
ggsave(file.path(paths$figures, "42_population_reliability_hist.png"), p_hist, width = 7, height = 5, dpi = 150)
message("Saved: 42_population_reliability_hist.png")

message("Stage 42 complete.")
