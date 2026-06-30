# Stage 11b: Posterior diagnostics for the stage-18 Student-t fit.
#
# The stage-18 post-sampling write_fit_outputs() call failed historically
# because the code looked for "mu" (Normal model parameter) in the t-model fit.
# That bug is resolved: summary_variables_for_model("18_real_lowrank_a_diag_b_t")
# correctly excludes mu (see src/stan_helpers.R:136). This script re-runs the
# diagnostic write and adds ESS/autocorrelation summaries by parameter group.
#
# Do NOT propose thinning as a runtime fix. Thinning reduces post-hoc draws
# but has zero effect on warmup/sampling cost (where >95% of runtime lives).

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/11b_stage18_diagnostics.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr", "posterior"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(posterior)
})

message("=== Stage 11b: Stage-18 t-model posterior diagnostics ===")

# ---- Load the stage-18 fit ------------------------------------------------------

fit <- load_cmdstan_fit(model_id = "18_real_lowrank_a_diag_b_t")

if (is.null(fit)) {
  stop("Stage-18 fit not found. Check fits/csv/18_real_lowrank_a_diag_b_t/")
}

# ---- Write standard diagnostic files (previously failed due to mu lookup) ------

message("Writing fit outputs for stage-18 t-model...")
write_fit_outputs(fit, "18_real_lowrank_a_diag_b_t")
message("Saved 18_real_lowrank_a_diag_b_t_diagnostic_summary.csv and sampler_diagnostics.csv")

# ---- Extended ESS / Rhat by parameter group ------------------------------------

message("Computing ESS / Rhat by parameter group (small parameters only)...")

small_vars <- c("nu", "lambda_a_diag", "Lambda_a", "psi_a",
                "sigma_b", "sigma_e", "prop_a", "prop_b", "prop_e")

ps <- fit$summary(small_vars, mean, sd, posterior::rhat,
                  posterior::ess_bulk, posterior::ess_tail)

# Assign parameter group for summarising
ps <- ps |>
  mutate(
    param_group = dplyr::case_when(
      variable == "nu"                          ~ "nu",
      grepl("^lambda_a_diag|^Lambda_a", variable) ~ "Lambda_a (loadings)",
      grepl("^psi_a",  variable)                ~ "psi_a (uniqueness)",
      grepl("^sigma_b", variable)               ~ "sigma_b (season scale)",
      grepl("^sigma_e", variable)               ~ "sigma_e (residual scale)",
      grepl("^prop_",  variable)                ~ "prop (variance fractions)",
      TRUE                                      ~ "other"
    )
  )

ess_by_group <- ps |>
  group_by(param_group) |>
  summarise(
    n_params    = n(),
    min_ess_bulk  = min(ess_bulk, na.rm = TRUE),
    med_ess_bulk  = median(ess_bulk, na.rm = TRUE),
    min_ess_tail  = min(ess_tail, na.rm = TRUE),
    med_ess_tail  = median(ess_tail, na.rm = TRUE),
    max_rhat      = max(rhat, na.rm = TRUE),
    n_rhat_gt_101 = sum(!is.na(rhat) & rhat > 1.01),
    .groups = "drop"
  ) |>
  arrange(min_ess_bulk)

write_csv(ess_by_group, file.path(paths$diagnostics, "18_posterior_health.csv"))
message("Saved 18_posterior_health.csv")

message("\n--- ESS by parameter group (stage-18 t-model) ---")
print(as.data.frame(ess_by_group))

# Compare to stage-10 Normal diagnostics
diag10_path <- file.path(paths$diagnostics, "10_real_lowrank_a_diag_b_diagnostic_summary.csv")
if (file.exists(diag10_path)) {
  diag10 <- read_csv(diag10_path, show_col_types = FALSE)
  message(sprintf("\nStage-10 (Normal): min_ess_bulk=%.1f, min_ess_tail=%.1f, max_rhat=%.4f",
                  diag10$min_ess_bulk, diag10$min_ess_tail, diag10$max_rhat))
  diag18_path <- file.path(paths$diagnostics, "18_real_lowrank_a_diag_b_t_diagnostic_summary.csv")
  if (file.exists(diag18_path)) {
    diag18 <- read_csv(diag18_path, show_col_types = FALSE)
    message(sprintf("Stage-18 (t):      min_ess_bulk=%.1f, min_ess_tail=%.1f, max_rhat=%.4f",
                    diag18$min_ess_bulk, diag18$min_ess_tail, diag18$max_rhat))
  }
}

# ---- Autocorrelation for worst-ESS parameters -----------------------------------

message("\nComputing autocorrelation for 5 worst-ESS parameters...")

worst5 <- ps |>
  filter(!is.na(ess_bulk)) |>
  arrange(ess_bulk) |>
  head(5) |>
  pull(variable)

message("5 worst ESS_bulk parameters: ", paste(worst5, collapse = ", "))

draws_worst <- fit$draws(variables = worst5, format = "draws_array")

acf_results <- lapply(worst5, function(v) {
  x <- as.vector(draws_worst[, , v])
  ac <- acf(x, lag.max = 20, plot = FALSE)
  data.frame(
    variable = v,
    lag      = ac$lag[, 1, 1],
    acf      = ac$acf[, 1, 1]
  )
})
acf_df <- do.call(rbind, acf_results)
write_csv(acf_df, file.path(paths$diagnostics, "18_worst_ess_autocorrelation.csv"))
message("Saved 18_worst_ess_autocorrelation.csv")

message("\n--- Stage-11b note on thinning ---")
message("ESS reported above reflects the effective sample size AFTER fitting.")
message("Thinning (keeping every kth draw) would reduce ESS further and has no")
message("effect on warmup or sampling cost. It is NOT a runtime solution.")
message("Runtime improvements must come from model restructuring (see Part B/E of rundown).")

message("\nStage 11b complete.")
