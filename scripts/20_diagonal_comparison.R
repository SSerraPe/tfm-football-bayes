# Stage 20: Compare diagonal Sigma_a model (stage 03) vs low-rank Sigma_a model (stage 10).
#
# Both models use the same data (real La Liga, N=4944, P=53).
# We compare them using:
#   1. LOO-CV (leave-one-player-season out) based on log_lik[n] already in generated quantities.
#   2. A side-by-side scatter of ICC_player (prop_a) per feature for both models.
#   3. A scatter of observed vs model-implied correlations for both models.
#
# Note on priors: stage 03 uses sigma_a ~ N(0, 0.7) (full player SD) while stage 10
# uses psi_a ~ N(0, 0.4) (uniqueness only) because the low-rank part already absorbs
# most player variance. The comparison is informative but not perfectly prior-matched;
# the LOO difference should be interpreted as an overall fit improvement, not purely
# the effect of adding factors.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/20_diagonal_comparison.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "loo", "ggplot2"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(loo)
})

theme_set(theme_bw(base_size = 11))

# ---- Load both fits ----------------------------------------------------------

fit_diag <- load_cmdstan_fit(
  fit_path = file.path(paths$fits, "03_real_diagonal_additive_fit.rds"),
  model_id = "03_real_diagonal_additive"
)

fit_lr <- load_cmdstan_fit(
  fit_path = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
  model_id = "10_real_lowrank_a_diag_b"
)

if (is.null(fit_diag)) stop("Stage 03 fit not found. Run stage 03 first.", call. = FALSE)
if (is.null(fit_lr))   stop("Stage 10 fit not found. Run stage 10 first.", call. = FALSE)

# ---- LOO-CV comparison -------------------------------------------------------

message("Computing LOO for diagonal model (stage 03) ...")
ll_diag <- fit_diag$draws("log_lik", format = "matrix")
loo_diag <- loo(ll_diag, cores = 1)

message("Computing LOO for low-rank model (stage 10) ...")
ll_lr <- fit_lr$draws("log_lik", format = "matrix")
loo_lr <- loo(ll_lr, cores = 1)

loo_comp <- loo_compare(loo_diag, loo_lr)
cat("\nLOO comparison (diagonal vs low-rank):\n")
print(loo_comp)

write.csv(as.data.frame(loo_comp), file.path(paths$tables, "20_loo_diagonal_vs_lowrank.csv"))
message("LOO comparison saved to outputs/tables/20_loo_diagonal_vs_lowrank.csv")

# ---- ICC_player scatter: diagonal vs low-rank --------------------------------
# Both models report prop_a[p] (player ICC per feature).
# We use the posterior summary csvs from stages 03 and 10.

icc_path_diag <- file.path(paths$tables, "03_real_diagonal_additive_posterior_summary.csv")
icc_path_lr   <- file.path(paths$tables, "10_real_lowrank_a_diag_b_posterior_summary.csv")

if (file.exists(icc_path_diag) && file.exists(icc_path_lr)) {
  sum_diag <- read_csv(icc_path_diag, show_col_types = FALSE)
  sum_lr   <- read_csv(icc_path_lr,   show_col_types = FALSE)

  # Extract prop_a per feature from both summaries
  extract_prop_a <- function(summary_df, P) {
    var_names <- sprintf("prop_a[%d]", seq_len(P))
    summary_df |> filter(variable %in% var_names) |>
      mutate(feature_idx = as.integer(sub(".*\\[([0-9]+)\\]", "\\1", variable))) |>
      select(feature_idx, mean) |>
      arrange(feature_idx)
  }

  model_objects <- readRDS(file.path(paths$processed, "model_objects.rds"))
  P <- length(model_objects$variable_names)

  prop_a_diag <- extract_prop_a(sum_diag, P)
  prop_a_lr   <- extract_prop_a(sum_lr,   P)

  icc_compare <- tibble(
    feature   = model_objects$variable_names[prop_a_diag$feature_idx],
    icc_diag  = prop_a_diag$mean,
    icc_lr    = prop_a_lr$mean
  )

  write_csv(icc_compare, file.path(paths$tables, "20_icc_diagonal_vs_lowrank.csv"))

  p_icc <- ggplot(icc_compare, aes(icc_diag, icc_lr)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.6) +
    geom_point(alpha = 0.7, size = 1.8) +
    ggrepel::geom_text_repel(
      data = icc_compare |> filter(abs(icc_lr - icc_diag) > 0.05),
      aes(label = feature), size = 2.5, max.overlaps = 20
    ) +
    labs(
      x     = "ICC_player — diagonal model (stage 03)",
      y     = "ICC_player — low-rank model (stage 10)",
      title = "Per-feature ICC_player: diagonal vs low-rank",
      subtitle = "Points above the diagonal: low-rank assigns more variance to the player effect"
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1))

  ggplot2_tryCatch <- tryCatch(
    ggsave(file.path(paths$figures, "20_icc_diagonal_vs_lowrank.png"),
           p_icc, width = 7, height = 6, dpi = 150),
    error = function(e) {
      # ggrepel may not be installed; fall back to plain geom_point
      p_icc2 <- ggplot(icc_compare, aes(icc_diag, icc_lr)) +
        geom_abline(slope = 1, intercept = 0, colour = "grey60") +
        geom_point(alpha = 0.7, size = 1.8) +
        labs(x = "ICC_player — diagonal (stage 03)", y = "ICC_player — low-rank (stage 10)",
             title = "ICC_player: diagonal vs low-rank")
      ggsave(file.path(paths$figures, "20_icc_diagonal_vs_lowrank.png"),
             p_icc2, width = 7, height = 6, dpi = 150)
    }
  )
  message("ICC scatter saved.")
}

# ---- Summary note ------------------------------------------------------------

best_model <- rownames(loo_comp)[1]
elpd_diff  <- round(loo_comp["model2", "elpd_diff"], 1)
se_diff    <- round(loo_comp["model2", "se_diff"], 1)

write_notes(
  c(
    "Stage 20: Diagonal vs Low-Rank model comparison",
    "",
    "Diagonal model: additive_diagonal.stan (stage 03)",
    "  Sigma_a = diag(sigma_a^2), Sigma_b = diag(sigma_b^2)",
    "",
    "Low-rank model: additive_lowrank_a_diag_b.stan (stage 10)",
    "  Sigma_a = Lambda Lambda' + diag(psi_a^2), Sigma_b = diag(sigma_b^2)",
    "",
    "LOO-CV results (log predictive density; higher is better):",
    capture.output(print(loo_comp)),
    "",
    paste0("Model preferred by LOO: ", best_model),
    paste0("ELPD difference (low-rank vs diagonal): ", elpd_diff, " ± ", se_diff),
    "",
    "See outputs/tables/20_loo_diagonal_vs_lowrank.csv and",
    "    outputs/figures/20_icc_diagonal_vs_lowrank.png"
  ),
  file.path(paths$notes, "20_diagonal_vs_lowrank_notes.txt")
)

message("Stage 20 complete.")
