# Stage 40 — Posterior predictive checks for the K=3 production fit (stage 28).
#
# Stage 19 ran the same two checks against the stale K=2/pre-rebuild fit (stage 10).
# This rewrites both against the current production model, per the standing
# constraint that everything in the Results chapter reflects the K=3 model
# described in Chapter 4 (methodology.tex), not any leftover K=2-era artifact.
#
#   1. Correlation PPC: does the model reproduce observed cross-feature correlations?
#      Model-implied off-diagonal covariance comes from Sigma_a = Lambda_a Lambda_a' +
#      diag(psi_a^2) (season/residual covariances are diagonal, contribute nothing
#      off-diagonal). Computed per posterior draw (mean + 90% CI), same method as
#      stage 19, but on Lambda_a/psi_a draws pulled from the stage-28 CSVs via the
#      column-pass extractor (avoids the as_cmdstan_fit() bottleneck on these ~2.6GB
#      x 4 chain files -- same reasoning as stages 29/31/38).
#   2. Residual distribution: e_{n,p} = Y - a_hat - b_hat vs Normal(0, sigma_e_hat),
#      using the K=3 A_mean/B_mean persisted by stage 29
#      (29_player_effect_means_k3.csv / 29_season_effect_means_k3.csv).
#
# Reads:  fits/csv/28_real_lowrank_a_diag_b_t_k3/.../*.csv
#         outputs/tables/29_player_effect_means_k3.csv, 29_season_effect_means_k3.csv
#         outputs/tables/28_real_lowrank_a_diag_b_t_k3_posterior_summary.csv
#         data/processed/model_objects.rds, Y_scaled.csv
# Writes: outputs/tables/40_correlation_ppc_k3.csv
#         outputs/figures/40_correlation_ppc_scatter_k3.png
#         outputs/figures/40_residual_distribution_check_k3.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/40_ppc_k3.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "tidyr"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(tidyr); library(ggplot2)
})

# Repointed to the corrected production fit (Execution step F): stage 34a, minutes-scaled,
# fixed phi=0.5, K*=3 confirmed unchanged by Execution step D. The posterior-mean sigma_e
# used below for the correlation-check denominator needs no further change: it already IS
# sigma_e,p, the reference-exposure residual scale (Definition def:phi), since the fit's
# posterior mean of the `sigma_e` parameter is exactly that quantity, not an average over
# observations at varying exposure.
MODEL_ID <- "34a_real_lowrank_a_diag_b_t_mv_k3"
K <- 3L
py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")

mo            <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names
P             <- length(feature_names)

Y_obs <- as.matrix(read_csv(file.path(paths$processed, "Y_scaled.csv"), show_col_types = FALSE))
colnames(Y_obs) <- feature_names

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
  if (ret != 0) stop("Python extraction failed for chain ", chain_id)
  read.csv(out_path, check.names = FALSE)
}
SCRATCHPAD <- file.path(tempdir(), "stage40_extract")
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

# ── Check 1: correlation PPC (Lambda_a + psi_a draws, small: 144 + 48 = 192 cols) ──

lambda_names <- paste0("Lambda_a.", rep(seq_len(P), K), ".", rep(seq_len(K), each = P))
psi_names    <- paste0("psi_a.", seq_len(P))
params       <- c(lambda_names, psi_names)

message("Discovering column indices for Lambda_a + psi_a (", length(params), " cols)...")
col_idx <- find_col_indices(csv_files[1], params)
missing <- names(col_idx)[is.na(col_idx)]
if (length(missing) > 0) stop("Not found in header: ", paste(head(missing, 5), collapse = ", "))

message("Extracting Lambda_a + psi_a draws from ", length(csv_files), " chains ...")
chain_dfs <- lapply(seq_along(csv_files), function(i) {
  out <- file.path(SCRATCHPAD, sprintf("chain%d.csv", i))
  message("  chain ", i, " ...")
  extract_batch(csv_files[i], i, col_idx, out)
})
draws_df <- do.call(rbind, chain_dfs)
n_draws  <- nrow(draws_df)
message("Extracted ", n_draws, " draws.")

lambda_mat_all <- as.matrix(draws_df[, lambda_names, drop = FALSE]); storage.mode(lambda_mat_all) <- "double"
psi_mat_all    <- as.matrix(draws_df[, psi_names, drop = FALSE]);    storage.mode(psi_mat_all)    <- "double"
rm(draws_df); gc()

# sigma_e and sigma_b posterior means (for total-variance denominator) -- from the
# existing small posterior summary CSV, no extraction needed.
ps <- read_csv(file.path(paths$tables, paste0(MODEL_ID, "_posterior_summary.csv")), show_col_types = FALSE)
get_mean_vec <- function(prefix) {
  ps |> filter(grepl(paste0("^", prefix, "\\["), variable)) |>
    mutate(p = as.integer(sub(paste0(prefix, "\\[(\\d+)\\]"), "\\1", variable))) |>
    arrange(p) |> pull(mean)
}
sigma_e_mean_vec <- get_mean_vec("sigma_e")
sigma_b_mean_vec <- get_mean_vec("sigma_b")
stopifnot(length(sigma_e_mean_vec) == P, length(sigma_b_mean_vec) == P)

message("Computing model-implied correlation matrix across ", n_draws, " draws ...")
implied_corr_draws <- array(NA_real_, dim = c(n_draws, P, P))
for (s in seq_len(n_draws)) {
  L_s   <- matrix(lambda_mat_all[s, ], P, K)
  psi_s <- psi_mat_all[s, ]
  Sig_a <- tcrossprod(L_s) + diag(psi_s^2)
  total_var <- diag(Sig_a) + sigma_b_mean_vec^2 + sigma_e_mean_vec^2
  denom <- sqrt(outer(total_var, total_var)); denom[denom < 1e-10] <- 1e-10
  implied_corr_draws[s, , ] <- Sig_a / denom
  if (s %% 1000 == 0) message("  draw ", s, "/", n_draws)
}

implied_corr_mean <- apply(implied_corr_draws, c(2, 3), mean)
implied_corr_lo    <- apply(implied_corr_draws, c(2, 3), quantile, 0.05)
implied_corr_hi    <- apply(implied_corr_draws, c(2, 3), quantile, 0.95)
rm(implied_corr_draws); gc()

obs_corr  <- cor(Y_obs)
idx_pairs <- which(lower.tri(obs_corr), arr.ind = TRUE)
plot_df_k3 <- tibble(
  model     = sprintf("K=%d (production)", K),
  feature_k = feature_names[idx_pairs[, 1]], feature_l = feature_names[idx_pairs[, 2]],
  obs_corr  = obs_corr[idx_pairs],
  pred_corr = implied_corr_mean[idx_pairs], pred_lo = implied_corr_lo[idx_pairs], pred_hi = implied_corr_hi[idx_pairs]
)

# K=0 comparison panel: under the diagonal player-covariance model (no shared factors),
# Sigma_a's off-diagonal entries are exactly zero for every posterior draw regardless of
# any parameter values -- Sigma_b and Sigma_e are diagonal too (Sec. sec:covariance), so
# every off-diagonal model-implied correlation is identically zero. This needs no fit and
# no draws: it is a closed-form consequence of the K=0 model structure.
plot_df_k0 <- tibble(
  model     = "K=0 (diagonal)",
  feature_k = feature_names[idx_pairs[, 1]], feature_l = feature_names[idx_pairs[, 2]],
  obs_corr  = obs_corr[idx_pairs],
  pred_corr = 0, pred_lo = 0, pred_hi = 0
)

plot_df <- bind_rows(plot_df_k0, plot_df_k3) |>
  mutate(model = factor(model, levels = c(sprintf("K=0 (diagonal)"), sprintf("K=%d (production)", K))))
write_csv(plot_df, file.path(paths$tables, "40_correlation_ppc_k3.csv"))
message(sprintf("Correlation PPC: r(obs, model-implied) = %.3f [K=%d]; K=0 model-implied is identically zero.",
  cor(plot_df_k3$obs_corr, plot_df_k3$pred_corr), K))

# Titles/subtitles are deliberately omitted (G2): all explanatory content -- what the
# panels are, what r is, why K=0's line is flat at zero -- lives in the LaTeX caption.
p_corr <- ggplot(plot_df, aes(obs_corr, pred_corr)) +
  geom_errorbar(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.15, width = 0) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, colour = "firebrick", linewidth = 0.7) +
  facet_wrap(~model, nrow = 1) +
  labs(x = "Observed sample correlation", y = "Model-implied correlation (posterior mean ± 90% CI)") +
  theme_minimal(base_size = 11)
ggsave(file.path(paths$figures, "40_correlation_ppc_scatter_k3.png"), p_corr, width = 10, height = 5, dpi = 150)
message("Saved: 40_correlation_ppc_scatter_k3.png")

# ── Check 2: residual distribution (uses K=3 A_mean/B_mean from stage 29) ─────

player_means_path <- file.path(paths$tables, "29_player_effect_means_k3.csv")
season_means_path <- file.path(paths$tables, "29_season_effect_means_k3.csv")
if (!file.exists(player_means_path) || !file.exists(season_means_path)) {
  stop("Run scripts/29_outlier_study.R first to generate the K=3 A_mean/B_mean tables.")
}
A_hat_df <- read_csv(player_means_path, show_col_types = FALSE)
B_hat_df <- read_csv(season_means_path, show_col_types = FALSE)

player_index <- mo$player_index
season_index <- mo$season_index

A_hat <- as.matrix(A_hat_df[, feature_names, drop = FALSE])
B_hat <- as.matrix(B_hat_df[, feature_names, drop = FALSE])

A_obs <- A_hat[player_index, , drop = FALSE]
B_obs <- B_hat[season_index, , drop = FALSE]
E_obs <- Y_obs - A_obs - B_obs

names(sigma_e_mean_vec) <- feature_names

icc_path <- file.path(paths$tables, "41_icc_summary_k3.csv")
if (file.exists(icc_path)) {
  icc_df <- read_csv(icc_path, show_col_types = FALSE)
  top_player <- head(icc_df$feature[order(-icc_df$icc_player)], 2)
  top_season <- head(icc_df$feature[order(-icc_df$icc_season)], 2)
  low_icc    <- head(icc_df$feature[order(icc_df$icc_player + icc_df$icc_season)], 2)
  selected_features <- unique(c(top_player, top_season, low_icc))
} else {
  selected_features <- feature_names[seq_len(6)]
}

resid_long <- as.data.frame(E_obs[, selected_features, drop = FALSE]) |>
  pivot_longer(everything(), names_to = "feature", values_to = "residual") |>
  mutate(sigma_e = sigma_e_mean_vec[feature], feature = factor(feature, levels = selected_features))

p_resid <- ggplot(resid_long, aes(residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "steelblue", alpha = 0.7) +
  geom_line(
    data = resid_long |> group_by(feature) |>
      reframe(x = seq(min(residual), max(residual), length.out = 200)) |>
      left_join(resid_long |> group_by(feature) |> summarise(se = unique(sigma_e)), by = "feature") |>
      mutate(y = dnorm(x, 0, se)),
    aes(x = x, y = y), colour = "firebrick", linewidth = 0.8
  ) +
  facet_wrap(~feature, scales = "free", ncol = 2) +
  # Title/subtitle deliberately omitted (G2): the LaTeX caption states what the red curve
  # is and that the fitted model is Student-t, not Normal.
  labs(x = "Residual  (Y - a_hat - b_hat)", y = "Density") +
  theme_minimal(base_size = 11)
ggsave(file.path(paths$figures, "40_residual_distribution_check_k3.png"), p_resid, width = 8, height = 10, dpi = 150)
message("Saved: 40_residual_distribution_check_k3.png")

message("Stage 40 complete.")
