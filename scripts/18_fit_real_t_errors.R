# Stage 18: Fit the low-rank Sigma_a + diagonal Sigma_b model with Student-t residuals.
#
# This is identical to stage 10 except:
#   - No global mean mu (data is already Z-scaled).
#   - Residuals are t(nu) instead of Normal.
#   - nu ~ Gamma(2, 0.1): a large posterior nu means normal errors are adequate;
#     a small nu (< ~10) signals meaningful heavy tails / outlier player-seasons.
#
# After fitting, we compare LOO-CV scores against the stage 10 normal-errors model.
# A positive elpd_diff favours the t-errors model.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/18_fit_real_t_errors.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr", "loo"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(loo)
})

# ---- Load the same Stan data as stage 10 ------------------------------------

model_objects_path <- file.path(paths$processed, "model_objects.rds")
Y_path             <- file.path(paths$processed, "Y_scaled.csv")
stan_data_path     <- file.path(paths$stan_data, "real_lowrank_a_diag_b_stan_data.rds")

if (!file.exists(stan_data_path)) stop("Missing stan data. Run stage 10 first.", call. = FALSE)

model_objects <- readRDS(model_objects_path)
stan_data     <- readRDS(stan_data_path)

# ---- Build PCA init for the t-errors model ----------------------------------
# build_pca_init() includes mu in the init list. The t-errors Stan model has no
# mu, so we remove it. We also seed nu at 30 (close to Normal) so the sampler
# starts in a sensible place.

base_init_fn <- build_pca_init(stan_data)

init_fn <- function(chain_id = 1) {
  init <- base_init_fn(chain_id)
  init$mu <- NULL    # t-errors model has no mu parameter
  init$nu <- 30.0    # start near Normal; sampler will pull toward the posterior
  init
}

# ---- Fit the t-errors model -------------------------------------------------

fit <- fit_cmdstan_model(
  stan_file  = file.path(paths$stan, "additive_lowrank_a_diag_b_t.stan"),
  stan_data  = stan_data,
  fit_path   = file.path(paths$fits, "18_real_lowrank_a_diag_b_t_fit.rds"),
  seed       = pipeline$seed + 180L,
  model_id   = "18_real_lowrank_a_diag_b_t",
  init       = init_fn
)

write_fit_outputs(fit, "18_real_lowrank_a_diag_b_t")

# ---- LOO comparison: normal errors (stage 10) vs t errors (stage 18) --------
# Both models computed log_lik[n] in generated quantities. LOO-CV estimates
# the expected log predictive density for new observations.

compare_loo <- function(fit_normal, fit_t) {
  if (is.null(fit_normal) || is.null(fit_t)) {
    message("One or both fits unavailable; skipping LOO comparison.")
    return(invisible(NULL))
  }

  ll_normal <- fit_normal$draws("log_lik", format = "matrix")
  ll_t      <- fit_t$draws("log_lik", format = "matrix")

  loo_normal <- loo(ll_normal, cores = 1)
  loo_t      <- loo(ll_t,      cores = 1)

  comparison <- loo_compare(loo_normal, loo_t)
  print(comparison)

  out_path <- file.path(paths$tables, "18_loo_normal_vs_t.csv")
  write.csv(as.data.frame(comparison), out_path)
  message("LOO comparison saved to ", out_path)

  list(loo_normal = loo_normal, loo_t = loo_t, comparison = comparison)
}

fit_normal <- load_cmdstan_fit(
  fit_path = file.path(paths$fits, "10_real_lowrank_a_diag_b_fit.rds"),
  model_id = "10_real_lowrank_a_diag_b"
)

loo_results <- compare_loo(fit_normal, fit)

# ---- Report posterior of nu -------------------------------------------------

if (!is.null(fit)) {
  nu_summary <- fit$summary("nu")
  cat("\nPosterior summary for nu (degrees of freedom):\n")
  print(nu_summary)

  write_notes(
    c(
      "Stage 18: Student-t residuals",
      "",
      "Stan model: stan/additive_lowrank_a_diag_b_t.stan",
      "Same player factor structure as stage 10 (low-rank Sigma_a + diagonal Sigma_b).",
      "Difference: residuals are Student-t(nu) instead of Normal.",
      "Prior on nu: Gamma(2, 0.1).",
      "",
      paste0("nu posterior mean:   ", round(nu_summary$mean, 2)),
      paste0("nu posterior median: ", round(nu_summary$median, 2)),
      paste0("nu 95% CI: [", round(nu_summary$q5, 2), ", ", round(nu_summary$q95, 2), "]"),
      "",
      "Interpretation:",
      "  nu >> 30: t errors nearly Normal; no gain over stage 10.",
      "  nu < 10:  meaningful heavy tails; outlier player-seasons present.",
      "",
      "LOO comparison saved to outputs/tables/18_loo_normal_vs_t.csv"
    ),
    file.path(paths$notes, "18_real_lowrank_a_diag_b_t_notes.txt")
  )
}

message("Stage 18 complete.")
