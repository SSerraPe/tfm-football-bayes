# Build an observed player-season panel from the raw Wyscout player-match CSV.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/00_build_longitudinal_dataset.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(required_base_packages)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

if (!file.exists(paths$raw)) {
  stop("Raw CSV not found at: ", paths$raw, call. = FALSE)
}

raw <- read_raw_player_match_csv(paths$raw)
schema <- schema_summary(raw)
panel <- build_player_season_panel(raw, pipeline$minimum_minutes_build)

season_counts <- panel |>
  count(season_id, name = "player_season_rows") |>
  arrange(season_id)

aggregation_summary <- tibble(
  item = c(
    "raw_rows",
    "player_season_rows",
    "players",
    "seasons",
    "competitions",
    "minimum_minutes_build",
    "eligible_player_season_rows"
  ),
  value = c(
    nrow(raw),
    nrow(panel),
    n_distinct(panel$player_id),
    n_distinct(panel$season_id),
    n_distinct(panel$competition_id),
    pipeline$minimum_minutes_build,
    sum(panel$is_min_minutes_eligible)
  )
)

aggregation_rules <- tibble(
  field_family = c("Identifiers", "Exposure", "Totals", "Per-90 rates", "Success rates", "Source averages/percentages"),
  rule = c(
    "Grouped by player_id, season_id, and competition_id; player_name is carried from the first non-missing value.",
    "observed_matches counts distinct match_id; minutes sums total_minutes_on_field.",
    "All numeric total_* columns are summed across observed lower-level rows.",
    "For each summed total_* count, per90_* = 90 * total / summed minutes.",
    "Selected rate_* variables are recomputed from summed successful/attempt count pairs.",
    "Original average_* and percent_* fields are retained as minute-weighted audit summaries, not as primary modelling variables."
  )
)

write_csv(panel, file.path(paths$interim, "longitudinal_player_seasons.csv"))
write_csv(panel, file.path(paths$processed, "longitudinal_player_seasons.csv"))
write_csv(schema, file.path(paths$tables, "00_raw_schema_summary.csv"))
write_csv(season_counts, file.path(paths$tables, "00_player_counts_by_season.csv"))
write_csv(aggregation_summary, file.path(paths$tables, "00_aggregation_summary.csv"))
write_csv(aggregation_rules, file.path(paths$tables, "00_aggregation_rules.csv"))

saveRDS(
  list(
    raw_schema = schema,
    season_counts = season_counts,
    aggregation_summary = aggregation_summary,
    aggregation_rules = aggregation_rules,
    minimum_minutes_build = pipeline$minimum_minutes_build
  ),
  file.path(paths$interim, "00_longitudinal_build_metadata.rds")
)

message("00_build_longitudinal_dataset complete: ", nrow(panel), " player-season rows.")
