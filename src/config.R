# Shared configuration for the additive Bayesian modelling pipeline.

required_base_packages <- c("readr", "dplyr", "tidyr", "tibble", "stringr")

paths <- list(
  raw = file.path(model_root, "data/raw/stg_wyscout__matches_players_stats_compLaLiga_season_14-15.csv"),
  interim = file.path(model_root, "data/interim"),
  processed = file.path(model_root, "data/processed"),
  stan_data = file.path(model_root, "stan_data"),
  fits = file.path(model_root, "fits"),
  tables = file.path(model_root, "outputs/tables"),
  diagnostics = file.path(model_root, "outputs/diagnostics"),
  figures = file.path(model_root, "outputs/figures"),
  notes = file.path(model_root, "outputs/notes"),
  scripts = file.path(model_root, "scripts"),
  stan = file.path(model_root, "stan"),
  notebooks = file.path(model_root, "notebooks")
)

invisible(lapply(paths[setdiff(names(paths), "raw")], dir.create, recursive = TRUE, showWarnings = FALSE))

pipeline <- list(
  seed = as.integer(Sys.getenv("BFA_SEED", unset = "20260513")),
  minimum_minutes_build = as.numeric(Sys.getenv("BFA_MIN_MINUTES_BUILD", unset = "300")),
  minimum_minutes_model = as.numeric(Sys.getenv("BFA_MIN_MINUTES", unset = "450")),
  missingness_threshold = as.numeric(Sys.getenv("BFA_MISSINGNESS_THRESHOLD", unset = "0.20")),
  sigma_floor = as.numeric(Sys.getenv("BFA_SIGMA_FLOOR", unset = "0.05"))
)

simulation <- list(
  rank_a = as.integer(Sys.getenv("BFA_SIM_RANK_A", unset = "2")),
  rank_b = as.integer(Sys.getenv("BFA_SIM_RANK_B", unset = "1")),
  seed = as.integer(Sys.getenv("BFA_SIM_SEED", unset = as.character(pipeline$seed + 400L)))
)

sampler <- list(
  run_stan = tolower(Sys.getenv("BFA_RUN_STAN", unset = "false")) %in% c("true", "1", "yes", "y"),
  reuse_fit = tolower(Sys.getenv("BFA_REUSE_FIT", unset = "true")) %in% c("true", "1", "yes", "y"),
  chains = as.integer(Sys.getenv("BFA_CHAINS", unset = "4")),
  parallel_chains = as.integer(Sys.getenv("BFA_PARALLEL_CHAINS", unset = "4")),
  iter_warmup = as.integer(Sys.getenv("BFA_ITER_WARMUP", unset = "1000")),
  iter_sampling = as.integer(Sys.getenv("BFA_ITER_SAMPLING", unset = "1000")),
  adapt_delta = as.numeric(Sys.getenv("BFA_ADAPT_DELTA", unset = "0.95")),
  max_treedepth = as.integer(Sys.getenv("BFA_MAX_TREEDEPTH", unset = "12")),
  refresh = as.integer(Sys.getenv("BFA_REFRESH", unset = "100"))
)

curated_feature_variables <- c(
  "per90_goals", "per90_assists", "per90_xg_shot", "per90_xg_assist",
  "per90_shots", "per90_head_shots", "rate_shots_on_target",
  "rate_goal_conversion", "per90_received_pass", "per90_passes",
  "per90_forward_passes", "per90_back_passes", "per90_lateral_passes",
  "per90_long_passes", "per90_smart_passes",
  "per90_progressive_passes", "per90_passes_to_final_third",
  "per90_key_passes", "per90_through_passes", "per90_crosses",
  "per90_shot_assists", "rate_successful_passes",
  "rate_successful_progressive_passes", "rate_successful_crosses",
  "per90_dribbles", "per90_progressive_run", "per90_accelerations",
  "rate_successful_dribbles", "per90_touch_in_box", "per90_losses",
  "per90_own_half_losses", "per90_dangerous_own_half_losses",
  "per90_recoveries", "per90_opponent_half_recoveries",
  "per90_dangerous_opponent_half_recoveries",
  "per90_interceptions", "per90_clearances",
  "per90_shots_blocked", "per90_defensive_duels", "per90_offensive_duels",
  "per90_aerial_duels", "per90_pressing_duels", "per90_loose_ball_duels",
  "per90_sliding_tackles", "per90_dribbles_against",
  "rate_defensive_duels_won", "rate_offensive_duels_won",
  "rate_aerial_duels_won", "rate_dribbles_against_won",
  "rate_successful_sliding_tackles", "per90_fouls",
  "per90_fouls_suffered", "per90_yellow_cards"
)
