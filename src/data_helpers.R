# Shared IO and longitudinal panel construction helpers.

read_raw_player_match_csv <- function(path) {
  data <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NULL"),
    progress = FALSE
  )

  keep_as_character <- c("player_name", "_last_modified_at", "_last_imported_at")
  for (col in setdiff(names(data), keep_as_character)) {
    values <- data[[col]]
    non_missing <- values[!is.na(values)]
    if (length(non_missing) == 0) next

    numeric_like <- all(stringr::str_detect(non_missing, "^-?[0-9]+(\\.[0-9]+)?$"))
    if (numeric_like) {
      data[[col]] <- readr::parse_double(values)
    }
  }

  data
}

schema_summary <- function(data) {
  tibble::tibble(
    column = names(data),
    class = vapply(data, function(x) paste(class(x), collapse = "/"), character(1)),
    n_missing = vapply(data, function(x) sum(is.na(x)), integer(1)),
    prop_missing = n_missing / nrow(data),
    n_distinct = vapply(data, function(x) dplyr::n_distinct(x, na.rm = TRUE), integer(1))
  )
}

rate_specs <- function() {
  tibble::tribble(
    ~rate_name, ~numerator, ~denominator,
    "rate_duels_won", "total_duels_won", "total_duels",
    "rate_defensive_duels_won", "total_defensive_duels_won", "total_defensive_duels",
    "rate_offensive_duels_won", "total_offensive_duels_won", "total_offensive_duels",
    "rate_aerial_duels_won", "total_aerial_duels_won", "total_aerial_duels",
    "rate_successful_passes", "total_successful_passes", "total_passes",
    "rate_successful_smart_passes", "total_successful_smart_passes", "total_smart_passes",
    "rate_successful_passes_to_final_third", "total_successful_passes_to_final_third", "total_passes_to_final_third",
    "rate_successful_crosses", "total_successful_crosses", "total_crosses",
    "rate_successful_forward_passes", "total_successful_forward_passes", "total_forward_passes",
    "rate_successful_back_passes", "total_successful_back_passes", "total_back_passes",
    "rate_successful_through_passes", "total_successful_through_passes", "total_through_passes",
    "rate_successful_key_passes", "total_successful_key_passes", "total_key_passes",
    "rate_successful_vertical_passes", "total_successful_vertical_passes", "total_vertical_passes",
    "rate_successful_long_passes", "total_successful_long_passes", "total_long_passes",
    "rate_successful_dribbles", "total_successful_dribbles", "total_dribbles",
    "rate_shots_on_target", "total_shots_on_target", "total_shots",
    "rate_goal_conversion", "total_goals", "total_shots",
    "rate_successful_progressive_passes", "total_successful_progressive_passes", "total_progressive_passes",
    "rate_successful_sliding_tackles", "total_successful_sliding_tackles", "total_sliding_tackles",
    "rate_dribbles_against_won", "total_dribbles_against_won", "total_dribbles_against",
    "rate_field_aerial_duels_won", "total_field_aerial_duels_won", "total_field_aerial_duels",
    "rate_gk_saves", "total_gk_saves", "total_gk_shots_against",
    "rate_gk_successful_exits", "total_gk_successful_exits", "total_gk_exits",
    "rate_gk_aerial_duels_won", "total_gk_aerial_duels_won", "total_gk_aerial_duels"
  )
}

build_player_season_panel <- function(raw, minimum_minutes) {
  required_columns <- c(
    "player_id", "match_id", "season_id", "competition_id", "total_minutes_on_field"
  )
  missing_required <- setdiff(required_columns, names(raw))
  if (length(missing_required) > 0) {
    stop("The raw file is missing required columns: ",
         paste(missing_required, collapse = ", "), call. = FALSE)
  }

  if (!"player_name" %in% names(raw)) {
    raw <- dplyr::mutate(raw, player_name = NA_character_)
  }

  numeric_columns <- names(raw)[vapply(raw, is.numeric, logical(1))]
  total_columns <- intersect(names(raw)[stringr::str_detect(names(raw), "^total_")], numeric_columns)
  average_columns <- intersect(names(raw)[stringr::str_detect(names(raw), "^average_")], numeric_columns)
  percent_columns <- intersect(names(raw)[stringr::str_detect(names(raw), "^percent_")], numeric_columns)

  player_seasons <- raw |>
    dplyr::group_by(player_id, season_id, competition_id) |>
    dplyr::summarise(
      player_name = first_non_missing_character(player_name),
      observed_matches = dplyr::n_distinct(match_id),
      raw_rows = dplyr::n(),
      dplyr::across(dplyr::all_of(total_columns), ~ sum(.x, na.rm = TRUE)),
      dplyr::across(
        dplyr::all_of(average_columns),
        ~ weighted_mean_minutes(.x, total_minutes_on_field),
        .names = "source_weighted_{.col}"
      ),
      dplyr::across(
        dplyr::all_of(percent_columns),
        ~ weighted_mean_minutes(.x, total_minutes_on_field),
        .names = "source_weighted_{.col}"
      ),
      .groups = "drop"
    )

  player_seasons <- dplyr::rename(player_seasons, minutes = total_minutes_on_field)
  if ("total_matches" %in% names(player_seasons)) {
    player_seasons <- dplyr::rename(player_seasons, source_total_matches = total_matches)
  }

  count_columns_for_per90 <- setdiff(
    names(player_seasons)[stringr::str_detect(names(player_seasons), "^total_")],
    c("total_minutes_tagged")
  )

  for (col in count_columns_for_per90) {
    player_seasons[[stringr::str_replace(col, "^total_", "per90_")]] <-
      safe_divide(player_seasons[[col]], player_seasons$minutes, 90)
  }

  specs <- rate_specs()
  for (row_id in seq_len(nrow(specs))) {
    spec <- specs[row_id, ]
    if (all(c(spec$numerator, spec$denominator) %in% names(player_seasons))) {
      player_seasons[[spec$rate_name]] <-
        safe_divide(player_seasons[[spec$numerator]], player_seasons[[spec$denominator]])
    }
  }

  player_seasons |>
    dplyr::mutate(
      minimum_minutes_threshold = minimum_minutes,
      is_min_minutes_eligible = minutes >= minimum_minutes
    ) |>
    dplyr::arrange(season_id, competition_id, player_id)
}
