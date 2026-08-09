# Stage 34 — Fit minutes-scaled Student-t model (K=3, fixed and estimated phi).
#
# Fits two variants of additive_lowrank_a_diag_b_t with observation-level residual scaling:
#   epsilon_{n,p} ~ t(nu, 0, sigma_e[p] * (m_ref / minutes_n)^phi)
#
#   Variant A (model_id "34a"): fixed phi = 0.5 (Poisson rate theory)
#              Stan file: additive_lowrank_a_diag_b_t_mv.stan
#   Variant B (model_id "34b"): estimated phi ~ N(0.5, 0.2)
#              Stan file: additive_lowrank_a_diag_b_t_mv_phi.stan
#
# Stage 30 diagnostic confirmed (p < 1e-17) that low-minute player-seasons are
# 2.4× over-represented in the worst residuals. This model addresses that directly.
#
# K is read from BFA_SIM_RANK_A (default 3, matching stage 28).
# Run both variants with:
#   BFA_RUN_STAN=true BFA_REUSE_FIT=false BFA_SIM_RANK_A=3 Rscript scripts/34_fit_minutes_scaled_t.R
#
# After fits: compare nu, share of low-minute worst residuals, and ELPD vs stage 28.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/34_fit_minutes_scaled_t.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble)
})

# K selection: explicit BFA_SIM_RANK_A takes priority; otherwise read K* from stage 32.
rank_a_env <- Sys.getenv("BFA_SIM_RANK_A", unset = "")
if (nchar(rank_a_env) > 0L) {
  rank_a <- as.integer(rank_a_env)
  message("Using K=", rank_a, " from BFA_SIM_RANK_A.")
} else {
  elpd_path <- file.path(paths$tables, "32_player_holdout_elpd.csv")
  if (file.exists(elpd_path)) {
    elpd_tbl <- readr::read_csv(elpd_path, show_col_types = FALSE)
    complete  <- elpd_tbl |> dplyr::filter(!is.na(elpd_mv))
    if (nrow(complete) > 0L) {
      rank_a <- complete$K[which.max(complete$elpd_mv)]
      message("Using K*=", rank_a, " from player-holdout ELPD (stage 32).")
    } else {
      rank_a <- 3L
      message("Stage 32 ELPD table has no results yet — defaulting to K=3.")
    }
  } else {
    rank_a <- 3L
    message("Stage 32 ELPD table not found — defaulting to K=3 (run stage 32 first).")
  }
}
message("Stage 34: minutes-scaled t-model at K=", rank_a)

# ── Load data ──────────────────────────────────────────────────────────────────

mo <- readRDS(file.path(paths$processed, "model_objects.rds"))
Y  <- mo$Y
minutes <- mo$metadata$minutes   # N-vector of playing minutes
m_ref   <- median(minutes)       # reference = sample median (1585 min)

message(sprintf("N=%d P=%d I=%d S=%d | minutes range %.0f–%.0f (median %.0f)",
  nrow(Y), ncol(Y), length(mo$player_levels), length(mo$season_levels),
  min(minutes), max(minutes), m_ref))

# ── Build stan_data (extends build_lowrank_a_diag_b_stan_data) ─────────────────

stan_data_base <- build_lowrank_a_diag_b_stan_data(mo, Y, rank_a)
stan_data_base$compute_log_lik <- 0L   # skip log_lik on initial fit; enable for LOO

stan_data_mv <- c(stan_data_base, list(
  minutes = as.numeric(minutes),
  m_ref   = m_ref
))

saveRDS(stan_data_mv,
  file.path(paths$stan_data,
            sprintf("real_lowrank_a_diag_b_t_mv_k%d_stan_data.rds", rank_a)))

# ── PCA init adapted for t-model (same as stage 28) ──────────────────────────

base_init_fn <- build_pca_init(stan_data_base)
init_fn_fixed <- function(chain_id = 1) {
  init    <- base_init_fn(chain_id)
  init$mu <- NULL
  init$nu <- 30.0   # start near Normal; let data pull toward heavier tails
  init
}
init_fn_phi <- function(chain_id = 1) {
  init     <- base_init_fn(chain_id)
  init$mu  <- NULL
  init$nu  <- 30.0
  init$phi <- 0.5   # start at Poisson-theory value
  init
}

# ── Fit variant A: fixed phi = 0.5 ────────────────────────────────────────────

fit_a_id   <- sprintf("34a_real_lowrank_a_diag_b_t_mv_k%d", rank_a)
fit_a_path <- file.path(paths$fits, paste0(fit_a_id, "_fit.rds"))
stan_file_a <- file.path(paths$stan, "additive_lowrank_a_diag_b_t_mv.stan")

fit_a <- fit_cmdstan_model(
  stan_file  = stan_file_a,
  stan_data  = stan_data_mv,
  fit_path   = fit_a_path,
  seed       = pipeline$seed + 340L,
  model_id   = fit_a_id,
  init       = init_fn_fixed
)
write_fit_outputs(fit_a, fit_a_id)

# ── Fit variant B: estimated phi ───────────────────────────────────────────────

fit_b_id   <- sprintf("34b_real_lowrank_a_diag_b_t_mv_phi_k%d", rank_a)
fit_b_path <- file.path(paths$fits, paste0(fit_b_id, "_fit.rds"))
stan_file_b <- file.path(paths$stan, "additive_lowrank_a_diag_b_t_mv_phi.stan")

fit_b <- fit_cmdstan_model(
  stan_file  = stan_file_b,
  stan_data  = stan_data_mv,
  fit_path   = fit_b_path,
  seed       = pipeline$seed + 341L,
  model_id   = fit_b_id,
  init       = init_fn_phi
)
write_fit_outputs(fit_b, fit_b_id)

# ── Summary: nu and phi posteriors ────────────────────────────────────────────

summarise_param <- function(fit, param, label) {
  if (is.null(fit)) { message(label, ": fit not available"); return(invisible(NULL)) }
  s <- tryCatch(fit$summary(param), error = function(e) NULL)
  if (!is.null(s)) {
    message(sprintf("%s: mean=%.3f  95%%CI=[%.3f, %.3f]",
                    label, s$mean, s$q5, s$q95))
  }
}

message("\n--- Stage 34 posterior summaries ---")
summarise_param(fit_a, "nu",  "Variant A (fixed phi=0.5)  nu")
summarise_param(fit_b, "nu",  "Variant B (estimated phi)  nu")
summarise_param(fit_b, "phi", "Variant B (estimated phi)  phi")

message("\nExpected outcome if minutes hypothesis is correct:")
message("  nu rises vs stage-28 (4.902), low-minute observations no longer dominate tails.")
message("  phi_hat close to 0.5 supports the Poisson rate interpretation.")
message("\nNext: re-run stage 30 diagnostic on new residuals to check low-minute over-representation.")
message("Stage 34 complete.")
