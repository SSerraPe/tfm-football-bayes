# Block C — LOO comparison: K=2 (stage 18) vs K=3 (stage 28) t-models.
#
# Uses generate_quantities(fitted_params = csv_paths) to avoid loading 9.6 GB
# of draws into R memory. The CmdStan binary streams the CSV draws, runs the
# generated quantities block with compute_log_lik=1, and writes compact output.
#
# Compare ΔELPD direction against Pathfinder K=2→3 gain of +12,125.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/C_block_c_loo_k2_vs_k3.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr", "loo"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(loo) })

stan_file <- file.path(paths$stan, "additive_lowrank_a_diag_b_t.stan")
check_packages("cmdstanr")
mod <- cmdstanr::cmdstan_model(stan_file)

GQ_DIR <- file.path(paths$fits, "csv", "block_c_gq")
dir.create(GQ_DIR, recursive = TRUE, showWarnings = FALSE)

run_gq_and_loo <- function(model_id, stan_data_path, gq_tag) {
  message("\n=== Running generate_quantities for ", model_id, " ===")

  stan_data <- readRDS(stan_data_path)
  stan_data$compute_log_lik <- 1L

  csv_files <- discover_cmdstan_csv_files(model_id)
  if (length(csv_files) == 0) stop("No CSV files for ", model_id)
  message("  Input CSVs: ", length(csv_files), " chains")
  message("  Stan data: N=", stan_data$N, " P=", stan_data$P, " Q_a=", stan_data$Q_a)

  gq_out_dir <- file.path(GQ_DIR, gq_tag)
  dir.create(gq_out_dir, recursive = TRUE, showWarnings = FALSE)

  message("  Running generate_quantities (passing CSV paths directly)...")
  t0 <- proc.time()[["elapsed"]]

  gq_fit <- mod$generate_quantities(
    fitted_params   = csv_files,
    data            = stan_data,
    seed            = pipeline$seed + 300L,
    parallel_chains = sampler$parallel_chains,
    output_dir      = gq_out_dir,
    output_basename = paste0(gq_tag, "_gq")
  )

  elapsed <- proc.time()[["elapsed"]] - t0
  message(sprintf("  generate_quantities done in %.0f s (%.1f min)", elapsed, elapsed / 60))

  # Extract log_lik: shape [draws, chains, N]
  message("  Extracting log_lik draws...")
  lik_draws <- gq_fit$draws("log_lik", format = "matrix")  # [draws_total, N]
  message("  log_lik dims: ", nrow(lik_draws), " draws × ", ncol(lik_draws), " obs")

  # PSIS-LOO
  message("  Computing PSIS-LOO...")
  r_eff <- loo::relative_eff(exp(lik_draws),
                              chain_id = rep(seq_len(sampler$chains),
                                             each = nrow(lik_draws) / sampler$chains))
  loo_res <- loo::loo(lik_draws, r_eff = r_eff)

  message("  ELPD = ", round(loo_res$estimates["elpd_loo", "Estimate"], 1),
          " ± ", round(loo_res$estimates["elpd_loo", "SE"], 1))

  pk <- loo_res$diagnostics$pareto_k
  message(sprintf("  Pareto-k: >0.5: %.1f%%  >0.7: %.1f%%",
                  100 * mean(pk > 0.5), 100 * mean(pk > 0.7)))

  list(loo = loo_res, n_draws = nrow(lik_draws), n_obs = ncol(lik_draws))
}

# ── Stage 28: K=3 ────────────────────────────────────────────────────────────

loo_k3 <- run_gq_and_loo(
  model_id       = "28_real_lowrank_a_diag_b_t_k3",
  stan_data_path = file.path(paths$stan_data, "real_lowrank_a_diag_b_k3_stan_data.rds"),
  gq_tag         = "k3"
)

# ── Stage 18: K=2 ────────────────────────────────────────────────────────────

loo_k2 <- run_gq_and_loo(
  model_id       = "18_real_lowrank_a_diag_b_t",
  stan_data_path = file.path(paths$stan_data, "real_lowrank_a_diag_b_stan_data.rds"),
  gq_tag         = "k2"
)

# ── Comparison ────────────────────────────────────────────────────────────────

message("\n=== LOO comparison: K=3 vs K=2 ===")
comp <- loo::loo_compare(loo_k2$loo, loo_k3$loo)
print(comp)

elpd_k2 <- loo_k2$loo$estimates["elpd_loo", "Estimate"]
elpd_k3 <- loo_k3$loo$estimates["elpd_loo", "Estimate"]
se_k2   <- loo_k2$loo$estimates["elpd_loo", "SE"]
se_k3   <- loo_k3$loo$estimates["elpd_loo", "SE"]

delta_elpd <- elpd_k3 - elpd_k2
se_delta   <- sqrt(se_k2^2 + se_k3^2)   # approximate

message(sprintf("  ELPD(K=2) = %.1f ± %.1f", elpd_k2, se_k2))
message(sprintf("  ELPD(K=3) = %.1f ± %.1f", elpd_k3, se_k3))
message(sprintf("  ΔELPD(K3 - K2) = %.1f ± %.1f (approx SE)", delta_elpd, se_delta))
message(sprintf("  Pathfinder K=2→3 gain: +12,125 (directional reference)"))
if (delta_elpd > 0) {
  message(sprintf("  Direction: K=3 BETTER (consistent with Pathfinder)"))
} else {
  message(sprintf("  Direction: K=2 BETTER or TIE (inconsistent with Pathfinder — flag in journal)"))
}

# ── Save ─────────────────────────────────────────────────────────────────────

result_tbl <- dplyr::tibble(
  model       = c("K2_stage18", "K3_stage28"),
  elpd_loo    = c(elpd_k2, elpd_k3),
  se_loo      = c(se_k2, se_k3),
  n_draws     = c(loo_k2$n_draws, loo_k3$n_draws),
  n_obs       = c(loo_k2$n_obs,   loo_k3$n_obs),
  pk_gt_0p5   = c(mean(loo_k2$loo$diagnostics$pareto_k > 0.5),
                  mean(loo_k3$loo$diagnostics$pareto_k > 0.5)),
  pk_gt_0p7   = c(mean(loo_k2$loo$diagnostics$pareto_k > 0.7),
                  mean(loo_k3$loo$diagnostics$pareto_k > 0.7)),
  delta_elpd  = c(NA_real_, delta_elpd),
  se_delta    = c(NA_real_, se_delta)
)
write_csv(result_tbl, file.path(paths$tables, "C_loo_k2_vs_k3.csv"))

message("\nSaved: C_loo_k2_vs_k3.csv")
message("Block C complete.")
