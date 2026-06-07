# Run one or more numbered pipeline stages in isolated Rscript processes.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else file.path(getwd(), "model", "run_pipeline.R")
model_root <- dirname(normalizePath(script_path, mustWork = TRUE))

stage_map <- c(
  "00" = "00_build_longitudinal_dataset.R",
  "01" = "01_read_clean_feature_engineering.R",
  "02" = "02_prepare_model_objects.R",
  "03" = "03_fit_real_diagonal_additive.R",
  "04" = "04_simulate_lowrank_additive_data.R",
  "05" = "05_fit_sim_diagonal_additive.R",
  "06" = "06_fit_sim_lowrank_recovery.R",
  "07" = "07_diagnostics_and_recovery_plots.R",
  "08" = "08_validate_all_models.R",
  "09" = "09_fit_sim_lowrank_a_diag_b.R",
  "10" = "10_fit_real_lowrank_a_diag_b.R",
  "11" = "11_loading_analysis_and_chains.R",
  "12" = "12_football_interpretation.R",
  "13" = "13_pca_postprocessing.R",
  "14" = "14_player_profiles.R",
  "15" = "15_season_profiles.R",
  "16" = "16_icc_analysis.R",
  "17" = "17_rank_selection_real.R",
  "18" = "18_fit_real_t_errors.R",
  "19" = "19_posterior_predictive_checks.R",
  "20" = "20_diagonal_comparison.R"
)

if (length(args) == 0L) {
  args <- names(stage_map)
}

unknown <- setdiff(args, names(stage_map))
if (length(unknown) > 0L) {
  stop("Unknown stage(s): ", paste(unknown, collapse = ", "),
       ". Valid stages are: ", paste(names(stage_map), collapse = ", "), call. = FALSE)
}

rscript <- file.path(R.home("bin"), "Rscript")
for (stage in args) {
  script <- file.path(model_root, "scripts", stage_map[[stage]])
  message("Running stage ", stage, ": ", basename(script))
  status <- system2(rscript, script)
  if (!identical(status, 0L)) {
    stop("Stage ", stage, " failed with exit status ", status, call. = FALSE)
  }
}

message("Pipeline complete: ", paste(args, collapse = ", "))
