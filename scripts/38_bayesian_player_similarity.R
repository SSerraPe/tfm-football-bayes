# Stage 38 — Fully Bayesian, per-posterior-draw player-similarity metric
# (production K=3 Student-t fit, stage 28).
#
# Existing similarity work (scripts/37_player_radar_similarity.R) uses a single
# point-estimate Euclidean distance on classical PCA scores. This script instead
# decomposes similarity into two posterior-draw-level quantities and reports a
# mean + 90% credible interval + cross-draw stability for each:
#
#   STYLE  — distance between rotated PCA factor scores f_i^(s) in the shared
#            K=3 factor space (same construction as scripts/13_pca_postprocessing.R's
#            G_pca_draws, but persisted here per-draw instead of only summarised).
#            f_i^(s) comes from the SHARED part of the player effect,
#            a_i = Lambda*eta_i + psi_a * z_u,i, i.e. it deliberately excludes the
#            idiosyncratic/uniqueness term.
#   LEVEL  — gap in modelled per90 output (A[i, per90_goals] and A[i, per90_xg_shot],
#            Z-scored/transformed units — A is the model's full player-effect
#            matrix, transformed parameter in stan/additive_lowrank_a_diag_b_t.stan)
#            between the same two players, per draw. This is the "how much do they
#            actually produce" axis, as opposed to "what shape is their profile."
#
# Two players can be STYLE-close (near-identical shared factor position) while
# being LEVEL-far apart (very different output) -- this is the formal, model-native
# answer to "can a player be stylistically similar to an outlier without matching
# their output": yes, and the model's own Sigma_a = Lambda*Lambda' + Psi_a split is
# exactly the mechanism that makes the distinction meaningful.
#
# NOTE: the dataset is La Liga only (2011/12-2022/23). Any "generational outlier"
# discussion (e.g. Haaland) is necessarily hypothetical/illustrative -- he is not
# in this data. Concrete case studies below use actual in-sample outliers instead
# (Messi, Cristiano Ronaldo, Neymar, etc.), already validated as case-study
# players in docs/professor/professor_summary.qmd.
#
# Does NOT call cmdstanr::as_cmdstan_fit() -- these CSVs are ~2.6GB x 4 chains and
# that call is a documented ~45 min / ~10GB RAM bottleneck (see CLAUDE.md,
# scripts/31_pca_with_ci.R, scripts/B_block_b_diagnostics_28.R). Uses the same
# Python column-pass pattern as B_block_b_diagnostics_28.R instead.
#
# Reads:  fits/csv/28_real_lowrank_a_diag_b_t_k3/.../*.csv (via discover_cmdstan_csv_files)
#         outputs/tables/29_player_factor_scores.csv (case-study player_index lookup)
#         outputs/tables/37_player_similarity_top10.csv (sanity-check baseline)
# Writes: outputs/tables/38_bayesian_similarity_case_studies.csv
#         outputs/tables/38_style_neighbor_stability.csv
#         outputs/figures/38_style_vs_level_scatter.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/38_bayesian_player_similarity.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(ggplot2)
})

MODEL_ID <- "28_real_lowrank_a_diag_b_t_k3"
K <- 3L
py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")

mo            <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names
P             <- length(feature_names)
I_players     <- length(mo$player_levels)
GOAL_IDX      <- match("per90_goals",   feature_names)  # 1
XG_IDX        <- match("per90_xg_shot", feature_names)  # 3
stopifnot(!is.na(GOAL_IDX), !is.na(XG_IDX))

csv_files <- discover_cmdstan_csv_files(MODEL_ID)
if (length(csv_files) == 0) stop("No CSV files found for ", MODEL_ID)
message("Found ", length(csv_files), " chain CSVs for ", MODEL_ID)

# ── Named case-study players (verified against outputs/tables/29_player_factor_scores.csv) ──
# NOTE: multiple "L. Suárez" rows exist (Wyscout id splits, e.g. transfer-related);
# player_index 431 (n_seasons = 8) is used here as the long-career Luis Suárez.
case_study_players <- tribble(
  ~player_index, ~player_name,
  63,  "L. Messi",
  39,  "Cristiano Ronaldo",
  583, "Neymar",
  33,  "Sergio Ramos",
  51,  "Gerard Piqué",
  431, "L. Suárez",
  38,  "K. Benzema"
)

# ── Column index discovery (single header parse) ─────────────────────────────

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
  if (ret != 0) stop("Python extraction failed for chain ", chain_id, " (", csv_path, ")")
  read.csv(out_path, check.names = FALSE)
}

SCRATCHPAD <- file.path(tempdir(), "stage38_extract")
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

# ── Pass 1: Lambda_a (P*K) + eta_a_raw (I*K) -- needed for style (all players) ──

lambda_names <- paste0("Lambda_a.",  rep(seq_len(P), K), ".", rep(seq_len(K), each = P))
eta_names    <- paste0("eta_a_raw.", rep(seq_len(I_players), K), ".", rep(seq_len(K), each = I_players))
pass1_params <- c(lambda_names, eta_names)

message("Pass 1: discovering column indices for Lambda_a + eta_a_raw (", length(pass1_params), " cols)...")
col_idx1 <- find_col_indices(csv_files[1], pass1_params)
missing1 <- names(col_idx1)[is.na(col_idx1)]
if (length(missing1) > 0) stop("Not found in header: ", paste(head(missing1, 5), collapse = ", "), " ...")

message("Pass 1: extracting from ", length(csv_files), " chains (this streams the full CSVs; expect a few minutes)...")
pass1_dfs <- lapply(seq_along(csv_files), function(i) {
  out <- file.path(SCRATCHPAD, sprintf("pass1_chain%d.csv", i))
  message("  chain ", i, " ...")
  extract_batch(csv_files[i], i, col_idx1, out)
})
pass1_df <- do.call(rbind, pass1_dfs)
n_iter   <- max(pass1_df[[".iteration"]])
n_chains <- max(pass1_df[[".chain"]])
n_draws  <- nrow(pass1_df)
message(sprintf("Pass 1 done: %d draws (%d iter x %d chains).", n_draws, n_iter, n_chains))

# Convert to plain matrices ONCE for fast row access in the per-draw loop below
# (repeated data.frame single-row indexing inside a ~4000-iteration loop is slow).
lambda_mat_all <- as.matrix(pass1_df[, lambda_names, drop = FALSE])
storage.mode(lambda_mat_all) <- "double"
eta_cols_by_q <- lapply(seq_len(K), function(q) paste0("eta_a_raw.", seq_len(I_players), ".", q))
eta_mat_all   <- as.matrix(pass1_df[, unlist(eta_cols_by_q), drop = FALSE])
storage.mode(eta_mat_all) <- "double"
rm(pass1_df); gc()

# ── Per-draw PCA post-processing (style scores), following scripts/13_pca_postprocessing.R ──

pca_from_B <- function(B, k) {
  eg <- eigen(B, symmetric = TRUE)
  vals <- pmax(eg$values[seq_len(k)], 0)
  list(U = eg$vectors[, seq_len(k), drop = FALSE], values = vals)
}

message("Computing reference orientation from posterior-mean Lambda_a ...")
B_ref <- matrix(0, P, P)
for (s in seq_len(n_draws)) {
  L_s <- matrix(lambda_mat_all[s, ], P, K)
  B_ref <- B_ref + tcrossprod(L_s) / n_draws
}
ref <- pca_from_B(B_ref, K)
L_ref <- ref$U %*% diag(sqrt(ref$values), K)
for (h in seq_len(K)) {
  jmax <- which.max(abs(L_ref[, h]))
  if (L_ref[jmax, h] < 0) L_ref[, h] <- -L_ref[, h]
}

message("Draw-wise PCA rotation of eta_a_raw -> style scores f_i^(s) (", n_draws, " draws x ", I_players, " players)...")
F_draws <- array(NA_real_, dim = c(n_draws, I_players, K))  # style scores per draw

for (s in seq_len(n_draws)) {
  L_s <- matrix(lambda_mat_all[s, ], P, K)
  Eta_s <- matrix(eta_mat_all[s, ], I_players, K)

  B_s <- tcrossprod(L_s)
  pca_s <- pca_from_B(B_s, K)
  safe_d <- pmax(pca_s$values, 1e-10)
  G_s <- (Eta_s %*% t(L_s)) %*% pca_s$U %*% diag(1 / sqrt(safe_d), K)

  for (h in seq_len(K)) {
    L_s_pca_h <- pca_s$U[, h] * sqrt(safe_d[h])
    if (sum(L_s_pca_h * L_ref[, h]) < 0) G_s[, h] <- -G_s[, h]
  }
  F_draws[s, , ] <- G_s
  if (s %% 500 == 0) message("  draw ", s, "/", n_draws)
}
saveRDS(F_draws, file.path(paths$fits, "38_style_score_draws.rds"))
rm(lambda_mat_all, eta_mat_all); gc()

# ── Style distance + stability, per case-study player vs all others ──────────

style_neighbors <- list()
style_summary_rows <- list()

for (r in seq_len(nrow(case_study_players))) {
  qi <- case_study_players$player_index[r]
  qn <- case_study_players$player_name[r]
  message("Style search for ", qn, " (player_index ", qi, ") ...")

  # Standardize each factor across players, per draw (so distance isn't dominated
  # by whichever factor happens to have larger raw scale in a given draw).
  d_style <- matrix(NA_real_, n_draws, I_players)
  for (s in seq_len(n_draws)) {
    Fs <- F_draws[s, , ]
    Fs_std <- scale(Fs)
    d_style[s, ] <- sqrt(rowSums(sweep(Fs_std, 2, Fs_std[qi, ], `-`)^2))
  }
  d_style[, qi] <- NA_real_  # exclude self

  mean_d <- colMeans(d_style, na.rm = TRUE)
  q05_d  <- apply(d_style, 2, quantile, probs = 0.05, na.rm = TRUE)
  q95_d  <- apply(d_style, 2, quantile, probs = 0.95, na.rm = TRUE)

  top10_idx <- order(mean_d)[seq_len(10)]
  # Stability: fraction of draws in which each of these top-10 is THE closest (rank 1)
  rank1_counts <- table(factor(apply(d_style, 1, which.min), levels = seq_len(I_players)))
  stability <- as.numeric(rank1_counts[top10_idx]) / n_draws

  style_neighbors[[qn]] <- tibble(
    query_player_index = qi, query_player_name = qn,
    neighbor_index = top10_idx,
    style_dist_mean = mean_d[top10_idx],
    style_dist_q05 = q05_d[top10_idx], style_dist_q95 = q95_d[top10_idx],
    p_closest_neighbor = stability
  )
}
style_neighbors_tbl <- bind_rows(style_neighbors)

player_names <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"), show_col_types = FALSE) |>
  select(player_index, player_id, player_name)
style_neighbors_tbl <- style_neighbors_tbl |>
  left_join(player_names, by = c("neighbor_index" = "player_index"))

write_csv(style_neighbors_tbl, file.path(paths$tables, "38_style_neighbor_stability.csv"))
message("Saved: 38_style_neighbor_stability.csv")

# ── Pass 2: A[i, goals] and A[i, xg_shot] for the small set of players that matter ──
# (case-study players + their top-10 style neighbors, ~70-80 unique players)

level_players <- sort(unique(c(case_study_players$player_index, style_neighbors_tbl$neighbor_index)))
message("Pass 2: extracting A[i, goals/xg_shot] for ", length(level_players), " players ...")

a_names <- c(
  paste0("A.", level_players, ".", GOAL_IDX),
  paste0("A.", level_players, ".", XG_IDX)
)
col_idx2 <- find_col_indices(csv_files[1], a_names)
missing2 <- names(col_idx2)[is.na(col_idx2)]
if (length(missing2) > 0) stop("Not found in header: ", paste(head(missing2, 5), collapse = ", "), " ...")

pass2_dfs <- lapply(seq_along(csv_files), function(i) {
  out <- file.path(SCRATCHPAD, sprintf("pass2_chain%d.csv", i))
  message("  chain ", i, " ...")
  extract_batch(csv_files[i], i, col_idx2, out)
})
pass2_df <- do.call(rbind, pass2_dfs)
stopifnot(nrow(pass2_df) == n_draws)
saveRDS(pass2_df, file.path(paths$fits, "38_pass2_A_goal_xg_draws.rds"))

goal_col <- function(i) paste0("A.", i, ".", GOAL_IDX)
xg_col   <- function(i) paste0("A.", i, ".", XG_IDX)

# ── Combine: style distance + level gap, per (query, neighbor) pair, per draw ──

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
      query_player_name   = qn,
      neighbor_player_name = nbrs$player_name[k],
      style_dist_mean = nbrs$style_dist_mean[k],
      style_dist_q05  = nbrs$style_dist_q05[k],
      style_dist_q95  = nbrs$style_dist_q95[k],
      p_closest_neighbor = nbrs$p_closest_neighbor[k],
      level_gap_mean = mean(level_gap),
      level_gap_q05  = quantile(level_gap, 0.05),
      level_gap_q95  = quantile(level_gap, 0.95)
    )
  }
}
case_studies_tbl <- bind_rows(case_rows) |> arrange(query_player_name, style_dist_mean)
write_csv(case_studies_tbl, file.path(paths$tables, "38_bayesian_similarity_case_studies.csv"))
message("Saved: 38_bayesian_similarity_case_studies.csv")

# ── Sanity check against existing point-estimate similarity table (stage 37) ──

classical_path <- file.path(paths$tables, "37_player_similarity_top10.csv")
if (file.exists(classical_path)) {
  classical <- read_csv(classical_path, show_col_types = FALSE)
  message("\n=== Sanity check vs classical (stage 37) top-1 neighbor ===")
  for (r in seq_len(nrow(case_study_players))) {
    qn <- case_study_players$player_name[r]
    new_top1 <- case_studies_tbl |> filter(query_player_name == qn) |> slice_min(style_dist_mean, n = 1)
    old_top1 <- classical |> filter(player_name == qn) |> slice_min(rank, n = 1)
    message(sprintf("  %s: new(Bayesian style) top-1 = %s | classical top-1 = %s",
                     qn,
                     if (nrow(new_top1) > 0) new_top1$neighbor_player_name[1] else "NA",
                     if (nrow(old_top1) > 0) old_top1$similar_player[1] else "NA"))
  }
}

# ── Figure: style distance vs level gap for the case-study pairs ─────────────

p <- ggplot(case_studies_tbl, aes(x = style_dist_mean, y = level_gap_mean, color = query_player_name)) +
  geom_errorbar(aes(ymin = level_gap_q05, ymax = level_gap_q95), width = 0, alpha = 0.3) +
  geom_errorbarh(aes(xmin = style_dist_q05, xmax = style_dist_q95), height = 0, alpha = 0.3) +
  geom_point(size = 2) +
  labs(
    title = "Style vs. level: top-10 style-neighbours of each case-study player",
    subtitle = "Bottom-left = near-identical style AND output. Bottom-right = same style, very different output level.",
    x = "Style distance (posterior mean, 90% CI)",
    y = "Level gap: |goals, xG per90| (posterior mean, 90% CI)",
    color = "Query player"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(paths$figures, "38_style_vs_level_scatter.png"), p, width = 9, height = 7, dpi = 150)
message("Saved: 38_style_vs_level_scatter.png")

message("\nStage 38 complete.")
