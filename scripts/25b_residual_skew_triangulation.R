# Task I — Triangulate Task A (heavy-tail features) vs Task C (skew candidates).
#
# Calibration: the expected p99 of |eps_std| under a correctly-specified
# t(nu=2.95) model is 3.37. Five of Task A's top-10 features by p99 fall
# at or below this level and are merely relative leaders among unremarkable
# residuals. Only features with observed p99 > 3.37 are genuinely underfitting.
#
# For each feature, classifies the cause of tail weight:
#   skew-driven     : skew > 1 AND in Task C's transform candidates
#   contamination   : excess p99 with keep-tier icc_player (real signal + outliers)
#   GK-contaminated : excess p99 on rate_defensive_duels_won and similar
#   unremarkable    : p99 below model expectation (only a relative leader)
#
# Output:
#   outputs/tables/25_residual_skew_triangulation.csv

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/25b_residual_skew_triangulation.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- Load inputs -------------------------------------------------------------

feat_summary <- read_csv(file.path(paths$tables, "22_residual_feature_summary.csv"),
                         show_col_types = FALSE)
skewness_df  <- read_csv(file.path(paths$tables, "23_feature_skewness.csv"),
                         show_col_types = FALSE)
icc_df       <- read_csv(file.path(paths$tables, "16_icc_summary.csv"),
                         show_col_types = FALSE) |>
  select(feature, icc_player, icc_season, icc_residual)

# Expected p99 of |eps_std| under t(nu=2.95): qt(0.995, df=2.95) / sqrt(nu/(nu-2))
nu_mean      <- 2.9508
scale_factor <- sqrt(nu_mean / (nu_mean - 2))
p99_expected <- qt(0.995, df = nu_mean) / scale_factor

message(sprintf("nu = %.4f  |  sqrt(nu/(nu-2)) = %.4f  |  Expected p99 = %.3f",
                nu_mean, scale_factor, p99_expected))

# ---- Join and classify -------------------------------------------------------

result <- feat_summary |>
  left_join(skewness_df, by = "feature") |>
  left_join(icc_df,      by = "feature") |>
  mutate(
    p99_excess = round(p99 - p99_expected, 3),
    genuinely_underfitting = p99 > p99_expected,
    verdict = dplyr::case_when(
      # Features NOT genuinely above model expectation at p99
      !genuinely_underfitting ~ "unremarkable_p99",
      # GK-driven: rate_defensive_duels_won (confirmed GK structural cluster)
      feature == "rate_defensive_duels_won" ~ "GK-contaminated",
      # Rate features that are drop candidates with confounded icc
      grepl("^rate_", feature) & icc_player < 0.20 ~ "drop_candidate_rate",
      # Keep-tier icc + excess p99 = contamination (not a transform target)
      icc_player >= 0.35 ~ "contamination-driven",
      # Season-driven: very high icc_season (league-wide tactical trend)
      !is.na(icc_season) & icc_season > 0.30 ~ "season-driven",
      # High skew + transform candidate
      !is.na(recommended_transform) & recommended_transform != "none" ~ "skew-driven",
      TRUE ~ "borderline"
    )
  ) |>
  select(feature, p99, p99_excess, max_abs_std, n_gt5sd, genuinely_underfitting,
         skewness, recommended_transform, icc_player, icc_season, verdict) |>
  arrange(desc(p99))

write_csv(result, file.path(paths$tables, "25_residual_skew_triangulation.csv"))
message("Saved 25_residual_skew_triangulation.csv")

# ---- Print calibration table -------------------------------------------------

message(sprintf("\n--- Calibration: expected p99 under t(%.4f) = %.3f ---", nu_mean, p99_expected))
message("Features with p99 genuinely above model expectation (p99_excess > 0):\n")
above <- result |> filter(genuinely_underfitting) |>
  select(feature, p99, p99_excess, max_abs_std, icc_player, verdict)
print(as.data.frame(above))

message("\nFeatures with p99 at or below model expectation (relative leaders only):\n")
below <- result |> filter(!genuinely_underfitting) |>
  select(feature, p99, p99_excess, icc_player, verdict)
print(as.data.frame(below))

# ---- Flag per90_smart_passes specifically ------------------------------------

if ("per90_smart_passes" %in% result$feature) {
  worst <- read_csv(file.path(paths$tables, "22_residual_worst_cells.csv"),
                    show_col_types = FALSE) |>
    filter(feature == "per90_smart_passes") |>
    arrange(desc(abs_std))
  message("\n--- per90_smart_passes: top cells (isolated vs. cluster?) ---")
  print(head(as.data.frame(worst |>
               select(player_name, season_id, eps_std)), 10))

  # Check if worst cells cluster by season
  worst_seasons <- worst |>
    count(season_id, name = "n_extreme_cells") |>
    arrange(desc(n_extreme_cells))
  message("Season cluster check (top 5 seasons by count of eps_std > 2):")
  extreme_by_season <- worst |>
    filter(abs_std > 2) |>
    count(season_id, sort = TRUE)
  print(head(as.data.frame(extreme_by_season), 5))
}

message("\nTask I complete.")
