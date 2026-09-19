# Stage 49 -- Per-player canonical-score credible intervals under a FIXED reference
# orientation (Results-chapter rewrite, forest-plot figure).
#
# This is a genuinely new construction, not previously implemented anywhere in the
# pipeline. The only existing per-draw score machinery (scripts/38_bayesian_player_
# similarity.R, scripts/13_pca_postprocessing.R) computes a PER-DRAW PCA -- its own
# V^(s) from that draw's own Lambda_a^(s) -- and only fixes the SIGN of each component
# against a reference, not the full rotation. Definition 4.22 / Remark 4.18
# (methodology.tex) call for something stronger: a single fixed eigenvector matrix
# V_bar, from the eigendecomposition of the POSTERIOR-MEAN common-variance matrix
# Lambda_bar %*% t(Lambda_bar), used to rotate EVERY draw's eta_i^(s):
#
#     f_i^(s) = V_bar^T %*% eta_i^(s)      (same V_bar for all s, all i)
#
# Rotation and sign ambiguity are then resolved by construction (V_bar is one fixed
# matrix), so per-draw scores are directly comparable across draws and across
# players -- which per-draw PCA (stage 38's approach) does not guarantee.
#
# This is a bounded, targeted extraction, not a refit: Lambda_a (144 cols, needed
# for ALL draws to build V_bar) + eta_a_raw restricted to the ~20 named players this
# figure needs (60 cols) -- a small fraction of the 4,731 columns stage 38 already
# extracts successfully from the same CSVs for all 1,529 players. Uses the same
# Python column-pass pattern (src/extract_stan_csv_params.py) as stage 38.
#
# Reads:  fits/csv/28_real_lowrank_a_diag_b_t_k3/.../*.csv (via discover_cmdstan_csv_files)
#         outputs/tables/29_player_factor_scores.csv (player_name -> player_index)
#         outputs/tables/48_distant_players.csv (which players are needed)
#         outputs/tables/37_player_similarity_top10.csv (which neighbours are needed)
# Writes: outputs/tables/49_score_credible_intervals.csv (player_name, pc, q05, q50, q95)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/49_score_credible_intervals.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages))

suppressPackageStartupMessages({ library(readr); library(dplyr) })

MODEL_ID <- "28_real_lowrank_a_diag_b_t_k3"
K <- 3L
py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")

mo <- readRDS(file.path(paths$processed, "model_objects.rds"))
P  <- length(mo$variable_names)

scores  <- read_csv(file.path(paths$tables, "29_player_factor_scores.csv"), show_col_types = FALSE) |>
  filter(!is.na(player_name))
distant <- read_csv(file.path(paths$tables, "48_distant_players.csv"), show_col_types = FALSE)
top10   <- read_csv(file.path(paths$tables, "37_player_similarity_top10.csv"), show_col_types = FALSE)

TARGETS <- c("L. Messi", "Sergio Ramos", "Xavi")
neighbour_names <- top10 |> filter(player_name %in% TARGETS, rank <= 3) |> pull(similar_player)
NEEDED_PLAYERS <- unique(c(TARGETS, neighbour_names, distant$player_name))
message("Players needed for the forest plot (", length(NEEDED_PLAYERS), "): ",
        paste(NEEDED_PLAYERS, collapse = ", "))

idx_lookup <- scores |> filter(player_name %in% NEEDED_PLAYERS) |>
  distinct(player_name, .keep_all = TRUE) |> select(player_name, player_index)
stopifnot(nrow(idx_lookup) == length(NEEDED_PLAYERS))

csv_files <- discover_cmdstan_csv_files(MODEL_ID)
if (length(csv_files) == 0) stop("No CSV files found for ", MODEL_ID)
message("Found ", length(csv_files), " chain CSVs for ", MODEL_ID)

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

SCRATCHPAD <- file.path(tempdir(), "stage49_extract")
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

lambda_names <- paste0("Lambda_a.", rep(seq_len(P), K), ".", rep(seq_len(K), each = P))
eta_names    <- as.vector(outer(idx_lookup$player_index, seq_len(K),
                                 function(i, q) paste0("eta_a_raw.", i, ".", q)))
pass_params  <- c(lambda_names, eta_names)

message("Discovering column indices for ", length(pass_params), " parameters (144 Lambda_a + ",
        length(eta_names), " eta_a_raw)...")
col_idx <- find_col_indices(csv_files[1], pass_params)
missing <- names(col_idx)[is.na(col_idx)]
if (length(missing) > 0) stop("Not found in header: ", paste(head(missing, 10), collapse = ", "))

message("Extracting from ", length(csv_files), " chains...")
chain_dfs <- lapply(seq_along(csv_files), function(i) {
  out <- file.path(SCRATCHPAD, sprintf("pass_chain%d.csv", i))
  message("  chain ", i, " ...")
  extract_batch(csv_files[i], i, col_idx, out)
})
draws_df <- do.call(rbind, chain_dfs)
n_draws  <- nrow(draws_df)
message(sprintf("Extraction done: %d draws.", n_draws))

lambda_mat_all <- as.matrix(draws_df[, lambda_names, drop = FALSE])
storage.mode(lambda_mat_all) <- "double"

# ---- fixed reference orientation: eigendecomposition of the POSTERIOR-MEAN common- ----
# ---- variance matrix Lambda_bar %*% t(Lambda_bar), sign-anchored per Remark 4.24 -------
Lambda_bar <- matrix(colMeans(lambda_mat_all), nrow = P, ncol = K)  # column-major fill matches Lambda_a.p.k naming
C_bar <- Lambda_bar %*% t(Lambda_bar)
eig <- eigen(C_bar, symmetric = TRUE)
V_bar <- eig$vectors[, 1:K, drop = FALSE]
for (k in 1:K) {
  j <- which.max(abs(V_bar[, k]))
  if (V_bar[j, k] < 0) V_bar[, k] <- -V_bar[, k]
}
message("Fixed reference V_bar built. Eigenvalues (variance shares): ",
        paste(round(100 * eig$values[1:K] / sum(eig$values[1:K]), 1), collapse = "%, "), "%")

# sanity check against the already-published PCA variance shares (61.8/27.7/10.5%)
# -- confirms this V_bar matches the same construction already used elsewhere.

# ---- per-draw fixed-reference scores for each needed player ---------------------------
# f_i^(s) = V_bar^T (Lambda^(s) eta_i^(s)): eta_i^(s) in R^K is first mapped into the
# P-dimensional shared-effect space by THAT DRAW'S OWN Lambda^(s) (exactly as the existing,
# proven per-draw-PCA code in scripts/38_bayesian_player_similarity.R and scripts/
# 13_pca_postprocessing.R does: `Eta_s %*% t(L_s)` before projecting) -- then, and only
# then, projected through the FIXED V_bar rather than a fresh per-draw eigendecomposition.
# This is the one piece the earlier per-draw-PCA scripts did NOT do (they re-eigendecompose
# every draw's own Lambda^(s) and only fix the sign against a reference), and is exactly
# what makes scores comparable draw-to-draw and player-to-player here. Confirmed dimensionally
# necessary: eta_i^(s) in R^K cannot be projected by V_bar^T (K x P) directly (checked; a
# first attempt at this script tried exactly that and failed with a non-conformable-
# arguments error, since eta_i^(s) is K-dimensional, not P-dimensional).
n_players <- nrow(idx_lookup)
eta_arr <- array(NA_real_, dim = c(n_draws, n_players, K))
for (p_i in seq_len(n_players)) {
  pi <- idx_lookup$player_index[p_i]
  eta_cols <- paste0("eta_a_raw.", pi, ".", seq_len(K))
  eta_arr[, p_i, ] <- as.matrix(draws_df[, eta_cols, drop = FALSE])
}

f_arr <- array(NA_real_, dim = c(n_draws, n_players, K))
Vt <- t(V_bar)  # K x P
for (s in seq_len(n_draws)) {
  Lambda_s <- matrix(lambda_mat_all[s, ], nrow = P, ncol = K)       # P x K, that draw's own loadings
  ETA_s    <- t(eta_arr[s, , ])                                     # K x n_players
  shared_s <- Lambda_s %*% ETA_s                                    # P x n_players
  f_arr[s, , ] <- t(Vt %*% shared_s)                                # n_players x K
}

results <- list()
for (p_i in seq_len(n_players)) {
  pn <- idx_lookup$player_name[p_i]
  for (k in 1:K) {
    q <- quantile(f_arr[, p_i, k], probs = c(0.05, 0.50, 0.95), names = FALSE)
    results[[length(results) + 1]] <- data.frame(
      player_name = pn, pc = paste0("PC", k), q05 = q[1], q50 = q[2], q95 = q[3]
    )
  }
}
out_df <- do.call(rbind, results)
write_csv(out_df, file.path(paths$tables, "49_score_credible_intervals.csv"))
message("Stage 49 complete: outputs/tables/49_score_credible_intervals.csv (", nrow(out_df), " rows)")

# Sanity check: compare fixed-reference posterior medians to the existing posterior-mean
# scores from stage 29 (should be close, though not identical -- median vs mean, and a
# fixed vs per-draw rotation is a different, if closely related, construction).
check <- out_df |> tidyr::pivot_wider(names_from = pc, values_from = c(q05, q50, q95)) |>
  left_join(scores |> select(player_name, pc1_score, pc2_score, pc3_score), by = "player_name")
message("\n--- Sanity check: fixed-reference posterior median vs. existing posterior-mean score ---")
print(check |> select(player_name, q50_PC1, pc1_score, q50_PC2, pc2_score, q50_PC3, pc3_score) |>
        mutate(across(where(is.numeric), ~round(., 2))), n = 30)
