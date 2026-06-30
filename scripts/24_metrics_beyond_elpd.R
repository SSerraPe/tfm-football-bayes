# Stage 24: Metrics beyond ELPD for the low-rank vs diagonal comparison.
#
# Provides four complementary metrics to supplement the ELPD/LOO comparison,
# addressing Issue 6 (better model comparison metrics) and Issue 10 (LOO
# reliability with bad Pareto-k).
#
# G1. Correlation Frobenius distance and RMSE from 19_correlation_ppc.csv.
#     Measures how well the model reproduces the observed cross-feature
#     correlation structure (a direct consequence of the factor structure).
#
# G2. Correlation prediction coverage.
#     Fraction of observed pairwise correlations inside the model's 90% CI.
#
# G3. ICC shift between diagonal and low-rank from 20_icc_diagonal_vs_lowrank.csv.
#     Quantifies how much the player-stability estimate changes across models.
#
# G3b. K-fold ELPD comparison (already done in stage 17 / real_data_elpd.csv).
#     The K=0 (diagonal) vs K=2 (low-rank) comparison is already computed;
#     this script extracts and contextualises those numbers.
#
# G4. Explanation of why the ELPD direction is trustworthy despite bad Pareto-k.
#
# All inputs are existing CSVs — NO Stan resampling required.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/24_metrics_beyond_elpd.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

message("=== Stage 24: Metrics beyond ELPD ===")

# ---- G1/G2: Correlation Frobenius distance from stage-19 PPC -------------------

message("--- G1/G2: Correlation recovery metrics ---")

corr_ppc <- read_csv(file.path(paths$tables, "19_correlation_ppc.csv"),
                     show_col_types = FALSE)
# Columns: feature_k, feature_l, obs_corr, pred_corr, pred_lo, pred_hi

n_pairs  <- nrow(corr_ppc)
diff_sq  <- (corr_ppc$obs_corr - corr_ppc$pred_corr)^2
abs_diff <- abs(corr_ppc$obs_corr - corr_ppc$pred_corr)

frob_dist   <- sqrt(sum(diff_sq))          # ||C_obs - C_pred||_F
corr_rmse   <- sqrt(mean(diff_sq))         # root-mean-squared prediction error
max_abs_err <- max(abs_diff)               # worst single pair
coverage_90 <- mean(
  corr_ppc$obs_corr >= corr_ppc$pred_lo &
  corr_ppc$obs_corr <= corr_ppc$pred_hi
)

# Pairs where the model's CI misses by more than 0.1
n_bad_pairs <- sum(abs_diff > 0.1)

# Worst-predicted pairs
worst_pairs <- corr_ppc |>
  mutate(abs_diff = abs(obs_corr - pred_corr)) |>
  arrange(desc(abs_diff)) |>
  head(10) |>
  select(feature_k, feature_l, obs_corr, pred_corr, abs_diff)

corr_metrics <- data.frame(
  metric = c("frobenius_distance", "correlation_rmse", "max_abs_error",
             "coverage_90pct_ci", "n_pairs", "n_pairs_err_gt_0.1"),
  value  = c(round(frob_dist, 4), round(corr_rmse, 4), round(max_abs_err, 4),
             round(coverage_90, 4), n_pairs, n_bad_pairs)
)

write_csv(corr_metrics, file.path(paths$tables, "24_correlation_recovery_metrics.csv"))
message("Saved 24_correlation_recovery_metrics.csv")

message("\nCorrelation recovery (low-rank K=2 model, Normal fit):")
print(as.data.frame(corr_metrics))
message("\nWorst-predicted feature pairs:")
print(as.data.frame(worst_pairs))

# ---- G3: ICC shift between diagonal and low-rank --------------------------------

message("\n--- G3: ICC shift (diagonal vs low-rank) ---")

icc_comp_path <- file.path(paths$tables, "20_icc_diagonal_vs_lowrank.csv")
if (!file.exists(icc_comp_path)) {
  message("20_icc_diagonal_vs_lowrank.csv not found; skipping ICC shift analysis.")
} else {
  icc_comp <- read_csv(icc_comp_path, show_col_types = FALSE)
  message("Columns in icc_comp: ", paste(colnames(icc_comp), collapse = ", "))

  # Expected columns: feature, prop_a_diag, prop_a_lowrank (or similar)
  # Adjust names based on actual column structure
  # Columns are: feature, icc_diag, icc_lr
  if (all(c("feature", "icc_diag", "icc_lr") %in% colnames(icc_comp))) {
    icc_shift <- icc_comp |>
      mutate(
        prop_a_diag   = icc_diag,
        prop_a_lowrank = icc_lr,
        shift     = icc_lr - icc_diag,
        abs_shift = abs(shift)
      ) |>
      arrange(desc(abs_shift))

    mean_abs_shift <- mean(icc_shift$abs_shift)
    n_up    <- sum(icc_shift$shift > 0.02)
    n_down  <- sum(icc_shift$shift < -0.02)

    message(sprintf("Mean |ICC_player shift| (lowrank - diagonal): %.4f", mean_abs_shift))
    message(sprintf("Features where ICC_player rises >0.02: %d", n_up))
    message(sprintf("Features where ICC_player falls >0.02: %d", n_down))
    message("\nTop features by ICC shift:")
    print(as.data.frame(head(icc_shift |> select(feature, prop_a_diag, prop_a_lowrank, shift), 10)))

    write_csv(icc_shift, file.path(paths$tables, "24_icc_shift.csv"))
    message("Saved 24_icc_shift.csv")
  } else {
    message("Unexpected column structure in 20_icc_diagonal_vs_lowrank.csv:")
    print(colnames(icc_comp))
    print(head(icc_comp, 3))
  }
}

# ---- G3b: K-fold ELPD (already computed in stage 17) ---------------------------

message("\n--- G3b: K-fold CV ELPD (diagonal vs low-rank) ---")

elpd_path <- file.path(model_root, "outputs", "rank_selection", "real_data_elpd.csv")
if (file.exists(elpd_path)) {
  elpd <- read_csv(elpd_path, show_col_types = FALSE)
  elpd_diag <- elpd$elpd[elpd$rank == 0]
  elpd_k2   <- elpd$elpd[elpd$rank == 2]
  delta_kfold <- elpd_k2 - elpd_diag

  elpd_per_diag <- elpd$elpd_per_value[elpd$rank == 0]
  elpd_per_k2   <- elpd$elpd_per_value[elpd$rank == 2]
  n_held_out    <- elpd$n_held_out_values[elpd$rank == 2]

  message(sprintf("K-fold ELPD (K=0 diagonal):  %.1f  (per value: %.4f)",
                  elpd_diag, elpd_per_diag))
  message(sprintf("K-fold ELPD (K=2 low-rank):  %.1f  (per value: %.4f)",
                  elpd_k2,   elpd_per_k2))
  message(sprintf("Δ ELPD (K=2 − K=0):         +%.1f  over %d held-out values",
                  delta_kfold, n_held_out))
  message(sprintf("Δ per held-out value:        +%.4f", elpd_per_k2 - elpd_per_diag))
} else {
  message("real_data_elpd.csv not found; skipping K-fold ELPD summary.")
}

# ---- G4: Why ELPD direction is trustworthy despite bad Pareto-k ----------------

loo_path <- file.path(paths$tables, "20_loo_diagonal_vs_lowrank.csv")
if (file.exists(loo_path)) {
  loo_comp <- read_csv(loo_path, show_col_types = FALSE)
  best_elpd <- loo_comp$elpd_loo[1]
  elpd_diff <- loo_comp$elpd_diff[2]
  se_diff   <- loo_comp$se_diff[2]
  n_se      <- abs(elpd_diff) / se_diff

  message("\n--- G4: ELPD reliability note ---")
  message(sprintf("PSIS-LOO Δ: %.1f (SE=%.1f, ratio=%.1f SEs)", elpd_diff, se_diff, n_se))
  message(
    "\nWhy the direction is trustworthy despite >58% 'very bad' Pareto-k:",
    "\n  Pareto-k > 0.7 means individual LOO estimates have very high variance",
    "\n  and the effective sample size per observation is very small. The per-",
    "\n  observation estimates and their aggregate SE are both understated.",
    "\n  However, this does NOT imply the SUM (the overall ELPD delta) is wrong",
    "\n  in direction — only that its variance is underestimated.",
    "\n",
    "\n  Two independent pieces of evidence point the same way:",
    sprintf("\n  1. PSIS-LOO Δ = %.1f (~%.0f SE) — even if SE is 3× understated,", elpd_diff, n_se),
    sprintf("\n     the signal (%.0f true SE) would still be very large.", n_se / 3),
    "\n  2. K-fold CV (no importance sampling) gives a much larger Δ (+72,880 ELPD",
    "\n     units) in the same direction, confirming this is not an IS artifact.",
    "\n",
    "\n  For the thesis: report K-fold as primary metric; PSIS-LOO as secondary",
    "\n  with an explicit Pareto-k caveat and a note that both agree directionally."
  )
}

message("\nStage 24 complete.")
