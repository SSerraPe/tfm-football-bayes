# Clean the longitudinal panel and construct curated non-GK per90/rate features.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/01_read_clean_feature_engineering.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

input_path <- file.path(paths$interim, "longitudinal_player_seasons.csv")
if (!file.exists(input_path)) {
  stop("Missing longitudinal dataset. Run model/scripts/00_build_longitudinal_dataset.R first.",
       call. = FALSE)
}

panel <- read_csv(input_path, show_col_types = FALSE, progress = FALSE, guess_max = Inf)
required_columns <- c("player_id", "season_id", "competition_id", "minutes", "observed_matches")
missing_required <- setdiff(required_columns, names(panel))
if (length(missing_required) > 0) {
  stop("The longitudinal panel is missing required columns: ",
       paste(missing_required, collapse = ", "), call. = FALSE)
}

if (!"player_name" %in% names(panel)) {
  panel <- panel |> mutate(player_name = NA_character_)
}

if ("total_gk_saves" %in% names(panel)) {
  n_pre_gk <- nrow(panel)
  panel <- panel |> filter(is.na(total_gk_saves) | total_gk_saves == 0)
  message("GK exclusion: removed ", n_pre_gk - nrow(panel), " GK player-seasons.")
} else {
  warning("Column 'total_gk_saves' not found in panel; GK players NOT excluded. ",
          "Run stage 00 to rebuild the interim file.")
}

eligible <- panel |>
  mutate(is_model_minutes_eligible = minutes >= pipeline$minimum_minutes_model) |>
  filter(is_model_minutes_eligible) |>
  arrange(season_id, competition_id, player_id)

metadata_columns <- intersect(
  c(
    "player_id", "player_name", "season_id", "competition_id",
    "observed_matches", "raw_rows", "minutes", "source_total_matches",
    "minimum_minutes_threshold", "is_min_minutes_eligible",
    "is_model_minutes_eligible"
  ),
  names(eligible)
)

metadata <- eligible |>
  select(all_of(metadata_columns)) |>
  mutate(row_id = row_number(), .before = 1)

candidate_feature_names <- names(eligible)[
  str_detect(names(eligible), "^(per90_|rate_)") &
    !str_detect(names(eligible), "gk_|goal_kicks")
]

candidate_features <- eligible |> select(all_of(candidate_feature_names))

missingness <- tibble(
  variable = names(candidate_features),
  n_missing = vapply(candidate_features, function(x) sum(is.na(x)), integer(1)),
  prop_missing = n_missing / nrow(candidate_features),
  n_observed = nrow(candidate_features) - n_missing
) |>
  arrange(desc(prop_missing), variable)

features_after_missingness <- missingness |>
  filter(prop_missing <= pipeline$missingness_threshold) |>
  pull(variable)

feature_matrix <- candidate_features |>
  select(all_of(features_after_missingness)) |>
  mutate(across(everything(), median_impute))

variance_summary <- tibble(
  variable = names(feature_matrix),
  variance = vapply(feature_matrix, stats::var, numeric(1)),
  n_distinct = vapply(feature_matrix, dplyr::n_distinct, integer(1))
)

features_after_variance <- variance_summary |>
  filter(is.finite(variance), variance > 1e-10, n_distinct > 2) |>
  pull(variable)

selected_variables <- curated_feature_variables[curated_feature_variables %in% features_after_variance]
if (length(selected_variables) < 10) {
  stop("The curated feature set left too few usable variables. Inspect model/src/config.R.",
       call. = FALSE)
}

final_features <- feature_matrix |> select(all_of(selected_variables))
if (any(!is.finite(as.matrix(final_features)))) {
  stop("Feature engineering produced non-finite model values.", call. = FALSE)
}

processed_feature_table <- bind_cols(metadata, final_features)

feature_audit <- tibble(variable = union(curated_feature_variables, candidate_feature_names)) |>
  mutate(
    is_candidate_per90_or_rate = variable %in% candidate_feature_names,
    in_curated_list = variable %in% curated_feature_variables,
    passed_missingness_filter = variable %in% features_after_missingness,
    passed_variance_filter = variable %in% features_after_variance,
    selected_for_model = variable %in% selected_variables,
    exclusion_reason = case_when(
      selected_for_model ~ "included",
      !is_candidate_per90_or_rate ~ "not a non-GK per90/rate candidate",
      !in_curated_list ~ "not in curated thesis feature set",
      !passed_missingness_filter ~ "missingness above threshold",
      !passed_variance_filter ~ "zero/near-zero variance or too few distinct values",
      TRUE ~ "not selected"
    )
  ) |>
  left_join(missingness, by = "variable") |>
  left_join(variance_summary, by = "variable") |>
  arrange(desc(selected_for_model), variable)

cleaning_summary <- tibble(
  item = c(
    "input_player_season_rows",
    "minimum_minutes_model",
    "eligible_player_season_rows",
    "candidate_features",
    "missingness_threshold",
    "features_after_missingness_filter",
    "features_after_variance_filter",
    "selected_features"
  ),
  value = c(
    nrow(panel),
    pipeline$minimum_minutes_model,
    nrow(eligible),
    length(candidate_feature_names),
    pipeline$missingness_threshold,
    length(features_after_missingness),
    length(features_after_variance),
    length(selected_variables)
  )
)

write_csv(eligible, file.path(paths$interim, "cleaned_player_seasons.csv"))
write_csv(metadata, file.path(paths$processed, "metadata_table.csv"))
write_csv(final_features, file.path(paths$processed, "model_feature_matrix_unscaled.csv"))
write_csv(processed_feature_table, file.path(paths$processed, "processed_feature_table.csv"))
write_csv(missingness, file.path(paths$tables, "01_feature_missingness_summary.csv"))
write_csv(variance_summary, file.path(paths$tables, "01_feature_variance_summary.csv"))
write_csv(feature_audit, file.path(paths$tables, "01_feature_inclusion_exclusion_audit.csv"))
write_csv(tibble(variable = selected_variables), file.path(paths$tables, "01_selected_variables.csv"))
write_csv(cleaning_summary, file.path(paths$tables, "01_cleaning_summary.csv"))

saveRDS(
  list(
    metadata_columns = metadata_columns,
    candidate_feature_names = candidate_feature_names,
    missingness_threshold = pipeline$missingness_threshold,
    curated_feature_variables = curated_feature_variables,
    selected_variables = selected_variables,
    minimum_minutes_model = pipeline$minimum_minutes_model
  ),
  file.path(paths$processed, "01_feature_engineering_metadata.rds")
)

message("01_read_clean_feature_engineering complete: ", length(selected_variables), " selected features.")
