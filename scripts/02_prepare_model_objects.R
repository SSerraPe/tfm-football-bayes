# Prepare model-ready matrices, player/season indices, maps, and Stan data.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/02_prepare_model_objects.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

feature_table_path <- file.path(paths$processed, "processed_feature_table.csv")
selected_variables_path <- file.path(paths$tables, "01_selected_variables.csv")
if (!file.exists(feature_table_path)) {
  stop("Missing processed feature table. Run stage 01 first.", call. = FALSE)
}
if (!file.exists(selected_variables_path)) {
  stop("Missing selected variables. Run stage 01 first.", call. = FALSE)
}

feature_table <- read_csv(feature_table_path, show_col_types = FALSE)
selected_variables <- read_csv(selected_variables_path, show_col_types = FALSE)$variable

required_metadata <- c("row_id", "player_id", "season_id", "competition_id", "minutes")
missing_metadata <- setdiff(required_metadata, names(feature_table))
if (length(missing_metadata) > 0) {
  stop("Processed feature table is missing metadata columns: ",
       paste(missing_metadata, collapse = ", "), call. = FALSE)
}

missing_features <- setdiff(selected_variables, names(feature_table))
if (length(missing_features) > 0) {
  stop("Processed feature table is missing selected feature columns: ",
       paste(missing_features, collapse = ", "), call. = FALSE)
}

metadata <- feature_table |>
  select(any_of(c(
    "row_id", "player_id", "player_name", "season_id", "competition_id",
    "observed_matches", "minutes"
  )))

Y_raw <- feature_table |>
  select(all_of(selected_variables)) |>
  as.matrix()

if (any(!is.finite(Y_raw))) {
  stop("The modelling matrix contains non-finite values after feature engineering.",
       call. = FALSE)
}

feature_means <- colMeans(Y_raw)
feature_sds <- apply(Y_raw, 2, stats::sd)
if (any(feature_sds <= 0 | !is.finite(feature_sds))) {
  bad <- names(feature_sds)[feature_sds <= 0 | !is.finite(feature_sds)]
  stop("Cannot scale zero-variance or non-finite features: ",
       paste(bad, collapse = ", "), call. = FALSE)
}

Y <- scale(Y_raw, center = feature_means, scale = feature_sds)
Y <- unclass(Y)
storage.mode(Y) <- "double"
colnames(Y) <- selected_variables

feature_matrix_audit <- audit_scaled_feature_matrix(Y, selected_variables)

player_levels <- sort(unique(metadata$player_id))
season_levels <- sort(unique(metadata$season_id))
player_index <- match(metadata$player_id, player_levels)
season_index <- match(metadata$season_id, season_levels)

row_mapping <- metadata |>
  mutate(
    model_row = row_number(),
    player_index = player_index,
    season_index = season_index,
    .before = 1
  )

player_map <- tibble(player_index = seq_along(player_levels), player_id = player_levels)
season_map <- tibble(season_index = seq_along(season_levels), season_id = season_levels)

player_recurrence <- row_mapping |>
  distinct(player_id, season_id, player_index, season_index) |>
  arrange(player_index, season_index) |>
  group_by(player_id, player_index) |>
  summarise(
    n_observed_seasons = n(),
    first_season_index = min(season_index),
    last_season_index = max(season_index),
    span_seasons = last_season_index - first_season_index + 1L,
    n_missing_within_span = span_seasons - n_observed_seasons,
    has_gap = n_missing_within_span > 0,
    gap_count = sum(diff(season_index) > 1L),
    .groups = "drop"
  ) |>
  mutate(recurrence_class = case_when(
    n_observed_seasons == 1L ~ "one-season",
    has_gap ~ "repeated-with-gap",
    TRUE ~ "consecutive-repeat"
  ))

season_counts <- row_mapping |>
  count(season_id, season_index, name = "player_season_rows") |>
  arrange(season_index)

scaling_parameters <- tibble(
  variable = selected_variables,
  center = as.numeric(feature_means),
  scale = as.numeric(feature_sds)
)

model_summary <- tibble(
  item = c(
    "N_player_seasons",
    "P_variables",
    "N_players",
    "N_seasons",
    "one_season_players",
    "repeated_players",
    "players_with_gaps",
    "feature_matrix_rank",
    "feature_matrix_max_abs_correlation"
  ),
  value = c(
    nrow(Y),
    ncol(Y),
    length(player_levels),
    length(season_levels),
    sum(player_recurrence$n_observed_seasons == 1L),
    sum(player_recurrence$n_observed_seasons > 1L),
    sum(player_recurrence$has_gap),
    feature_matrix_audit$summary$empirical_rank,
    feature_matrix_audit$summary$max_abs_correlation
  )
)

model_objects <- list(
  Y = Y,
  Y_raw = Y_raw,
  metadata = metadata,
  row_mapping = row_mapping,
  player_index = as.integer(player_index),
  season_index = as.integer(season_index),
  player_levels = player_levels,
  season_levels = season_levels,
  player_map = player_map,
  season_map = season_map,
  variable_names = selected_variables,
  scaling_parameters = scaling_parameters,
  player_recurrence = player_recurrence,
  season_counts = season_counts,
  feature_matrix_audit = feature_matrix_audit,
  model_summary = model_summary
)

stan_data <- build_diagonal_stan_data(model_objects)

write_csv(as_tibble(Y, .name_repair = ~ selected_variables), file.path(paths$processed, "Y_scaled.csv"))
write_csv(row_mapping, file.path(paths$processed, "model_row_mapping.csv"))
write_csv(player_map, file.path(paths$processed, "player_map.csv"))
write_csv(season_map, file.path(paths$processed, "season_map.csv"))
write_csv(scaling_parameters, file.path(paths$processed, "scaling_parameters.csv"))
write_csv(player_recurrence, file.path(paths$tables, "02_player_recurrence_diagnostics.csv"))
write_csv(season_counts, file.path(paths$tables, "02_player_counts_by_season.csv"))
write_csv(feature_matrix_audit$summary, file.path(paths$tables, "02_feature_matrix_audit_summary.csv"))
write_csv(feature_matrix_audit$high_correlation_pairs, file.path(paths$tables, "02_feature_high_correlation_pairs.csv"))
write_csv(model_summary, file.path(paths$tables, "02_model_object_summary.csv"))
saveRDS(model_objects, file.path(paths$processed, "model_objects.rds"))
saveRDS(stan_data, file.path(paths$stan_data, "real_diagonal_additive_stan_data.rds"))

write_notes(
  c(
    "Irregular Panel Treatment",
    "",
    "The modelling unit is an observed player-season row.",
    "Missing player-seasons are not imputed or represented as missing responses.",
    "Player and season indices only map rows that are observed in the Wyscout-derived panel."
  ),
  file.path(paths$notes, "02_irregular_panel_note.txt")
)

message("02_prepare_model_objects complete: N=", nrow(Y), ", P=", ncol(Y),
        ", players=", length(player_levels), ", seasons=", length(season_levels), ".")
