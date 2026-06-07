# Stage 19: Posterior predictive checks for the stage 10 model.
#
# Two complementary checks:
#
# 1. Correlation PPC: does the model reproduce the observed cross-variable
#    correlations? For each posterior draw we compute the model-implied
#    correlation from Sigma_a = Lambda_a Lambda_a' + diag(psi_a^2). The
#    off-diagonal entries of Sigma_a give the model-implied cross-feature
#    covariance; dividing by total marginal SDs gives the implied correlation.
#    We compare to the sample correlation matrix of the observed Y.
#
# 2. Residual distribution check: compute e_{n,p} = Y - a_hat - b_hat and
#    compare the empirical residual distribution to Normal(0, sigma_e_hat).
#    If the model is well-specified, the residuals should look Normal with
#    the right scale.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/19_posterior_predictive_checks.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "tidyr", "posterior"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

theme_set(theme_bw(base_size = 11))

# ---- Load data and model objects --------------------------------------------

model_objects <- readRDS(file.path(paths$processed, "model_objects.rds"))
Y_path        <- file.path(paths$processed, "Y_scaled.csv")
Y_obs         <- as.matrix(read_csv(Y_path, show_col_types = FALSE))
colnames(Y_obs) <- model_objects$variable_names

player_index  <- model_objects$player_index
season_index  <- model_objects$season_index
feature_names <- model_objects$variable_names
P             <- ncol(Y_obs)
N             <- nrow(Y_obs)
Q_a           <- simulation$rank_a

# Posterior mean player and season effects saved by stages 14 and 15
player_means_path <- file.path(paths$tables, "14_player_effect_means.csv")
season_means_path <- file.path(paths$tables, "15_season_effect_means.csv")

if (!file.exists(player_means_path) || !file.exists(season_means_path)) {
  stop("Run stages 14 and 15 first to generate player and season effect means.", call. = FALSE)
}

player_means <- read_csv(player_means_path, show_col_types = FALSE)
season_means <- read_csv(season_means_path, show_col_types = FALSE)

# ---- Check 1: Cross-variable correlation PPC --------------------------------
# Model-implied off-diagonal covariance between features k and l comes
# entirely from the player covariance Sigma_a = Lambda_a Lambda_a' + Psi_a.
# (Season and residual covariances are diagonal.)
#
# We use Lambda_a and psi_a draws from the stage 10 posterior.

fit_10 <- load_cmdstan_fit(
  fit_path = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
  model_id = "10_real_lowrank_a_diag_b"
)

if (!is.null(fit_10)) {
  message("Extracting Lambda_a and psi_a draws ...")
  n_draws_use <- 500L

  lambda_mat <- posterior::as_draws_matrix(fit_10$draws("Lambda_a"))
  psi_mat    <- posterior::as_draws_matrix(fit_10$draws("psi_a"))
  sigma_e_mat <- posterior::as_draws_matrix(fit_10$draws("sigma_e"))
  sigma_b_mat <- posterior::as_draws_matrix(fit_10$draws("sigma_b"))

  set.seed(42L)
  use_idx <- sample(seq_len(nrow(lambda_mat)), min(n_draws_use, nrow(lambda_mat)))

  Lambda_var_names <- as.vector(
    outer(seq_len(P), seq_len(Q_a), function(p, q) sprintf("Lambda_a[%d,%d]", p, q))
  )
  psi_var_names   <- sprintf("psi_a[%d]",   seq_len(P))
  sigma_e_names   <- sprintf("sigma_e[%d]", seq_len(P))
  sigma_b_names   <- sprintf("sigma_b[%d]", seq_len(P))

  # For each draw: Sigma_a = L L' + diag(psi^2); total Var(y_k) = Sigma_a[k,k] + sigma_b[k]^2 + sigma_e[k]^2
  # Model-implied correlation corr_model(k,l) = Sigma_a[k,l] / sqrt(Var(y_k) * Var(y_l))

  implied_corr_draws <- array(NA_real_, dim = c(length(use_idx), P, P))

  message("Computing model-implied correlation matrix across ", length(use_idx), " draws ...")
  for (idx in seq_along(use_idx)) {
    s     <- use_idx[idx]
    L_s   <- matrix(as.numeric(lambda_mat[s, Lambda_var_names]), P, Q_a)
    psi_s <- as.numeric(psi_mat[s, psi_var_names])
    se_s  <- as.numeric(sigma_e_mat[s, sigma_e_names])
    sb_s  <- as.numeric(sigma_b_mat[s, sigma_b_names])

    Sig_a <- L_s %*% t(L_s) + diag(psi_s^2)
    total_var <- diag(Sig_a) + sb_s^2 + se_s^2

    denom <- sqrt(outer(total_var, total_var))
    denom[denom < 1e-10] <- 1e-10
    implied_corr_draws[idx, , ] <- Sig_a / denom
  }

  implied_corr_mean <- apply(implied_corr_draws, c(2, 3), mean)
  implied_corr_lo   <- apply(implied_corr_draws, c(2, 3), quantile, 0.05)
  implied_corr_hi   <- apply(implied_corr_draws, c(2, 3), quantile, 0.95)

  # Observed sample correlation of Y
  obs_corr <- cor(Y_obs)

  # Scatter plot: observed vs model-implied (off-diagonal only)
  idx_pairs     <- which(lower.tri(obs_corr), arr.ind = TRUE)
  plot_df <- tibble(
    feature_k  = feature_names[idx_pairs[, 1]],
    feature_l  = feature_names[idx_pairs[, 2]],
    obs_corr   = obs_corr[idx_pairs],
    pred_corr  = implied_corr_mean[idx_pairs],
    pred_lo    = implied_corr_lo[idx_pairs],
    pred_hi    = implied_corr_hi[idx_pairs]
  )

  write_csv(plot_df, file.path(paths$tables, "19_correlation_ppc.csv"))

  p_corr <- ggplot(plot_df, aes(obs_corr, pred_corr)) +
    geom_errorbar(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.15, width = 0) +
    geom_point(alpha = 0.35, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, colour = "firebrick", linewidth = 0.7) +
    labs(
      x     = "Observed sample correlation",
      y     = "Model-implied correlation (posterior mean ± 90% CI)",
      title = "Cross-variable correlation: observed vs model-implied",
      subtitle = "Points near the diagonal indicate the model captures the correlation structure"
    )

  ggsave(file.path(paths$figures, "19_correlation_ppc_scatter.png"),
         p_corr, width = 7, height = 6, dpi = 150)
  message("Correlation PPC plot saved.")
}

# ---- Check 2: Residual distribution -----------------------------------------
# Residuals e_{n,p} = Y[n,p] - a_hat[player_index[n], p] - b_hat[season_index[n], p]
# should look like Normal(0, sigma_e_hat[p]) if the model is adequate.

message("Building residual matrix ...")

# Player means: wide format (player_id, feature1, feature2, ...)
player_id_col   <- names(player_means)[1]
common_features <- intersect(feature_names, names(player_means))

if (length(common_features) < P) {
  message("Warning: player effect means available for ", length(common_features),
          " of ", P, " features.")
}

A_hat <- as.matrix(player_means[, common_features, drop = FALSE])
rownames(A_hat) <- as.character(player_means[[player_id_col]])

# Season means: long format (season, feature, b_mean, ...)
# Pivot to wide: one row per season, one column per feature
season_wide <- season_means |>
  filter(feature %in% common_features) |>
  select(season, feature, b_mean) |>
  pivot_wider(names_from = feature, values_from = b_mean)
season_id_col <- "season"
B_hat <- as.matrix(season_wide[, common_features, drop = FALSE])
rownames(B_hat) <- as.character(season_wide[[season_id_col]])

player_ids_obs <- as.character(model_objects$player_levels[player_index])
season_ids_obs <- as.character(model_objects$season_levels[season_index])

# Align rows to observation order
A_obs <- A_hat[player_ids_obs, common_features]
B_obs <- B_hat[season_ids_obs, common_features]
Y_sub <- Y_obs[, common_features]
E_obs <- Y_sub - A_obs - B_obs

# Posterior mean sigma_e per feature
if (!is.null(fit_10)) {
  sigma_e_mean <- colMeans(sigma_e_mat[, sigma_e_names])
  names(sigma_e_mean) <- feature_names
} else {
  sigma_e_mean <- rep(1, P)
  names(sigma_e_mean) <- feature_names
}

# Plot residuals for a selection of features (cover different ICC levels)
icc_path <- file.path(paths$tables, "16_icc_summary.csv")
if (file.exists(icc_path)) {
  icc_df <- read_csv(icc_path, show_col_types = FALSE)
  # Pick 6 features: top-2 ICC_player, top-2 ICC_season, 2 low-ICC
  top_player  <- head(icc_df$feature[order(-icc_df$icc_player)], 2)
  top_season  <- head(icc_df$feature[order(-icc_df$icc_season)], 2)
  low_icc     <- head(icc_df$feature[order(icc_df$icc_player + icc_df$icc_season)], 2)
  selected_features <- unique(c(top_player, top_season, low_icc))
} else {
  selected_features <- common_features[seq(1, min(6, length(common_features)))]
}

selected_features <- intersect(selected_features, common_features)

resid_long <- as.data.frame(E_obs[, selected_features]) |>
  pivot_longer(everything(), names_to = "feature", values_to = "residual") |>
  mutate(
    sigma_e = sigma_e_mean[feature],
    feature = factor(feature, levels = selected_features)
  )

p_resid_clean <- ggplot(resid_long, aes(residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "steelblue", alpha = 0.7) +
  geom_line(
    data = resid_long |>
      group_by(feature) |>
      reframe(x = seq(min(residual), max(residual), length.out = 200)) |>
      left_join(resid_long |> group_by(feature) |> summarise(se = unique(sigma_e)), by = "feature") |>
      mutate(y = dnorm(x, 0, se)),
    aes(x = x, y = y), colour = "firebrick", linewidth = 0.8
  ) +
  facet_wrap(~feature, scales = "free", ncol = 2) +
  labs(
    x = "Residual  (Y − â_i − b̂_j)",
    y = "Density",
    title = "Residual distributions vs Normal(0, σ_e)",
    subtitle = "Red curve: posterior mean Normal(0, σ_e). Histogram: observed residuals."
  )

ggsave(file.path(paths$figures, "19_residual_distribution_check.png"),
       p_resid_clean, width = 8, height = 10, dpi = 150)
message("Residual distribution plot saved.")

write_notes(
  c(
    "Stage 19: Posterior predictive checks",
    "",
    "Check 1: Cross-variable correlation",
    "  Model-implied correlation from Sigma_a = Lambda_a Lambda_a' + diag(psi_a^2).",
    "  Compare to observed sample correlation of Y (already Z-scaled).",
    "  Result: see outputs/figures/19_correlation_ppc_scatter.png",
    "  Table:  outputs/tables/19_correlation_ppc.csv",
    "",
    "Check 2: Residual distributions",
    "  Residuals = Y - posterior_mean(a_i) - posterior_mean(b_j).",
    "  Compared to Normal(0, posterior_mean(sigma_e)) for selected features.",
    "  Result: see outputs/figures/19_residual_distribution_check.png",
    "",
    paste0("Selected features for residual check: ",
           paste(selected_features, collapse = ", "))
  ),
  file.path(paths$notes, "19_posterior_predictive_checks_notes.txt")
)

message("Stage 19 complete.")
