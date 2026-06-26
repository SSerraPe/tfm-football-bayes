profile = Sys.getenv("PLAYER_MODEL_SUMMARY_PROFILE", "thorough")
results_dir = file.path("..", "results", "reproduction")

main = readRDS(file.path(results_dir, paste0("main_simulation_", profile, ".rds")))
supplementary = readRDS(
  file.path(results_dir, paste0("supplementary_comparisons_", profile, ".rds"))
)

simulation_elpd = aggregate(
  elpd_per_value ~ true_rank + rank,
  main$cv_summary,
  mean
)
names(simulation_elpd)[3] = "mean_elpd_per_reconstructed_value"

supplementary_elpd = supplementary$summary[
  ,
  c("true_rank", "task", "rank", "elpd")
]

write.csv(
  main$runs,
  file.path(results_dir, paste0("note_simulation_study_runs_", profile, ".csv")),
  row.names = FALSE
)
write.csv(
  simulation_elpd,
  file.path(results_dir, paste0("note_simulation_study_elpd_", profile, ".csv")),
  row.names = FALSE
)
write.csv(
  supplementary_elpd,
  file.path(results_dir, paste0("note_supplementary_comparisons_elpd_", profile, ".csv")),
  row.names = FALSE
)

cat("Wrote the three note tables for profile:", profile, "\n")
