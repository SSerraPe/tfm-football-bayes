# Block F — Marginalized model sim recovery with I=1000 (5× original I=200).
#
# Unit test: verify S_i=2 off-diagonal P×P block of Cov_i = Sigma_a exactly.
# If PASS: run NUTS with I=1000 and report coverage improvement over I=200 (78%).
# True parameters identical to stage 27 for comparability.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/B_block_f_marginalized_i1000.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr", "posterior", "MASS"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(posterior); library(MASS)
})

SIM_P  <- 10L; SIM_K  <- 2L; SIM_I  <- 1000L; SIM_S  <- 10L
SIM_SEED <- pipeline$seed + 270L   # same seed as stage 27 for reproducibility

# ── 1. True parameters (identical to stage 27) ──────────────────────────────

set.seed(SIM_SEED)

L_raw <- matrix(0, SIM_P, SIM_K)
L_raw[c(1,3,5,7), 1] <- c( 1.2,  0.8, -0.5,  0.4)
L_raw[c(2,4,6),   1] <- c( 0.3, -0.6,  0.7)
L_raw[c(2,4,6,8,9,10), 2] <- c(0.9, -0.4, 0.5, 0.6, -0.3, 0.7)
Lambda_a_true <- make_plt(L_raw)

psi_a_true   <- c(0.5, 0.4, 0.6, 0.3, 0.5, 0.7, 0.4, 0.5, 0.6, 0.3)
sigma_b_true <- rep(0.4, SIM_P)
sigma_e_true <- rep(0.6, SIM_P)
Sigma_a_true <- tcrossprod(Lambda_a_true) + diag(psi_a_true^2)

# ── 2. Unit test: S_i=2 off-diagonal block = Sigma_a ────────────────────────

D_e_true <- diag(sigma_e_true^2)
Cov_2P <- rbind(
  cbind(Sigma_a_true + D_e_true, Sigma_a_true),
  cbind(Sigma_a_true, Sigma_a_true + D_e_true)
)
offdiag_block <- Cov_2P[1:SIM_P, (SIM_P + 1L):(2L * SIM_P)]
if (!isTRUE(all.equal(offdiag_block, Sigma_a_true, tolerance = 1e-12))) {
  stop("UNIT TEST FAILED: S_i=2 off-diagonal P×P block != Sigma_a")
}
message("UNIT TEST PASS: S_i=2 off-diagonal P×P block of Cov_i equals Sigma_a exactly.")
message("Proceeding to I=1000 sim recovery...")

# ── 3. Simulate data (I=1000) ────────────────────────────────────────────────

set.seed(SIM_SEED + 1L)

S_count  <- sample(1L:5L, SIM_I, replace = TRUE)
N_total  <- sum(S_count)

player_idx <- rep(seq_len(SIM_I), times = S_count)
season_idx <- unlist(lapply(S_count, function(si) sort(sample(SIM_S, si, replace = FALSE))))

z_b_true <- matrix(rnorm(SIM_S * SIM_P), SIM_S, SIM_P)
z_b_true <- sweep(z_b_true, 2, colMeans(z_b_true), "-")
B_true   <- sweep(z_b_true, 2, sigma_b_true, "*")

Y_sim <- matrix(NA_real_, N_total, SIM_P)
row_r <- 1L
for (i in seq_len(SIM_I)) {
  a_i <- drop(mvrnorm(1L, rep(0, SIM_P), Sigma_a_true))
  for (ss in seq_len(S_count[i])) {
    s <- season_idx[row_r]
    Y_sim[row_r, ] <- a_i + B_true[s, ] + rnorm(SIM_P, 0, sigma_e_true)
    row_r <- row_r + 1L
  }
}

message(sprintf("Simulated: N=%d  I=%d  S=%d  P=%d  K=%d", N_total, SIM_I, SIM_S, SIM_P, SIM_K))
message(sprintf("  Seasons per player: min=%d max=%d mean=%.1f",
                min(S_count), max(S_count), mean(S_count)))

# ── 4. Stan data ─────────────────────────────────────────────────────────────

n_free <- SIM_K * SIM_P - SIM_K * (SIM_K + 1L) %/% 2L

stan_data_sim <- list(
  N               = N_total,
  P               = SIM_P,
  I               = SIM_I,
  S               = SIM_S,
  Q_a             = SIM_K,
  N_lambda_a_free = n_free,
  Y               = unname(Y_sim),
  player_index    = as.integer(player_idx),
  season_index    = as.integer(season_idx),
  S_count         = as.integer(S_count),
  sigma_floor     = pipeline$sigma_floor
)

# ── 5. Init ──────────────────────────────────────────────────────────────────

lambda_a_diag0 <- diag(Lambda_a_true[1:SIM_K, 1:SIM_K, drop = FALSE])
lambda_a_free0 <- numeric(n_free)
pos <- 1L
for (h in seq_len(SIM_K)) {
  for (j in (h + 1L):SIM_P) {
    lambda_a_free0[pos] <- Lambda_a_true[j, h]
    pos <- pos + 1L
  }
}

init_fn <- function(chain_id = 1) {
  list(
    lambda_a_diag = abs(lambda_a_diag0 + rnorm(SIM_K, 0, 0.05)) + 0.05,
    lambda_a_free = lambda_a_free0 + rnorm(n_free, 0, 0.05),
    psi_a   = psi_a_true   + abs(rnorm(SIM_P, 0, 0.05)),
    z_b     = z_b_true + matrix(rnorm(SIM_S * SIM_P, 0, 0.1), SIM_S, SIM_P),
    sigma_b = sigma_b_true + abs(rnorm(SIM_P, 0, 0.05)),
    sigma_e = sigma_e_true + abs(rnorm(SIM_P, 0, 0.05))
  )
}

# ── 6. Fit ───────────────────────────────────────────────────────────────────

stan_file  <- file.path(paths$stan, "player_season_lowrank_marginalized.stan")
model_id   <- "27b_sim_lowrank_marginalized_i1000"
output_dir <- file.path(paths$fits, "csv", model_id, format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mod <- cmdstanr::cmdstan_model(stan_file)

message("\nFitting marginalized model I=1000 (warmup=500, sampling=500) ...")
fit <- mod$sample(
  data            = stan_data_sim,
  seed            = SIM_SEED,
  chains          = 4L,
  parallel_chains = 4L,
  iter_warmup     = 500L,
  iter_sampling   = 500L,
  adapt_delta     = 0.90,
  max_treedepth   = 10L,
  refresh         = 100L,
  output_dir      = output_dir,
  output_basename = model_id,
  init            = init_fn
)

# ── 7. Diagnostics ────────────────────────────────────────────────────────────

diag_sum  <- fit$diagnostic_summary(quiet = TRUE)
la_sum    <- fit$summary("Lambda_a")
message("\n=== Sampler diagnostics ===")
message("  Max Rhat (Lambda_a): ", round(max(la_sum$rhat, na.rm = TRUE), 3))
message("  Min ESS  (Lambda_a): ", round(min(la_sum$ess_bulk, na.rm = TRUE), 0))
message("  Divergences total:   ", sum(diag_sum$num_divergent))

# ── 8. Coverage ───────────────────────────────────────────────────────────────

post_sum <- fit$summary(c("Lambda_a", "psi_a", "sigma_b", "sigma_e"))

Lambda_a_rows <- expand.grid(p = 1:SIM_P, q = 1:SIM_K)
Lambda_a_true_vec <- Lambda_a_true[cbind(Lambda_a_rows$p, Lambda_a_rows$q)]

true_vals <- tibble::tibble(
  variable = c(
    paste0("Lambda_a[", Lambda_a_rows$p, ",", Lambda_a_rows$q, "]"),
    paste0("psi_a[",   1:SIM_P, "]"),
    paste0("sigma_b[", 1:SIM_P, "]"),
    paste0("sigma_e[", 1:SIM_P, "]")
  ),
  true_value = c(Lambda_a_true_vec, psi_a_true, sigma_b_true, sigma_e_true)
)

coverage_tbl <- post_sum |>
  dplyr::select(variable, q5, q95, mean, rhat, ess_bulk) |>
  dplyr::left_join(true_vals, by = "variable") |>
  dplyr::mutate(
    covered_90 = !is.na(true_value) & true_value >= q5 & true_value <= q95,
    abs_bias   = abs(mean - true_value)
  )

n_testable    <- sum(!is.na(coverage_tbl$true_value))
n_covered     <- sum(coverage_tbl$covered_90, na.rm = TRUE)
coverage_rate <- n_covered / n_testable

message("\n=== Coverage check (90% posterior CI) ===")
message("  Parameters tested: ", n_testable)
message(sprintf("  Covered: %d / %d = %.1f%% (I=200 baseline: 78.0%%)",
                n_covered, n_testable, 100 * coverage_rate))

misses <- coverage_tbl |>
  dplyr::filter(!is.na(true_value), !covered_90) |>
  dplyr::arrange(desc(abs_bias))

if (nrow(misses) > 0) {
  message("  Misses (", nrow(misses), "):")
  for (i in seq_len(min(nrow(misses), 10))) {
    m <- misses[i, ]
    message(sprintf("    %s: true=%.3f  CI=[%.3f, %.3f]  bias=%.3f",
                    m$variable, m$true_value, m$q5, m$q95, m$abs_bias))
  }
}

verdict <- if (coverage_rate >= 0.85) "PASS" else "FAIL"
message("\n  VERDICT: ", verdict, " (coverage ", if (coverage_rate >= 0.85) ">=" else "<", " 85%)")

# ── 9. Save ───────────────────────────────────────────────────────────────────

write_csv(coverage_tbl, file.path(paths$tables, "27b_marginalized_i1000_coverage.csv"))

note_lines <- c(
  "# Simulation Recovery: Marginalized Model — I=1000 (Block F)",
  "",
  sprintf("P=%d  K=%d  I=%d  S=%d  N=%d", SIM_P, SIM_K, SIM_I, SIM_S, N_total),
  sprintf("Seasons per player: 1-5 (mean %.1f)", mean(S_count)),
  "",
  "## Unit test",
  "S_i=2 off-diagonal P×P block of Cov_i equals Sigma_a = ΛΛ' + Ψ_a: PASS",
  "",
  "## Coverage (90% CI)",
  sprintf("Covered: %d / %d = %.1f%%", n_covered, n_testable, 100 * coverage_rate),
  "Baseline (I=200, stage 27): 78.0%",
  sprintf("Verdict: %s", verdict),
  "",
  "## Sampler",
  sprintf("Max Rhat (Lambda_a): %.3f", max(la_sum$rhat, na.rm = TRUE)),
  sprintf("Min ESS  (Lambda_a): %.0f", min(la_sum$ess_bulk, na.rm = TRUE)),
  sprintf("Divergences total:   %d",   sum(diag_sum$num_divergent)),
  "",
  sprintf("_Written by scripts/B_block_f_marginalized_i1000.R on %s_", Sys.Date())
)

writeLines(note_lines, file.path(paths$notes, "note_marginalized_i1000_recovery.md"))
message("\nSaved: 27b_marginalized_i1000_coverage.csv, note_marginalized_i1000_recovery.md")
message("Block F complete.")
