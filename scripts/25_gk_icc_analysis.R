# Task I2 — GK player impact across all 53 features.
#
# Confirms that goalkeeper players are present in the model data (not excluded —
# only GK-specific feature columns were filtered in stage 01). Flags GK rows via
# total_gk_saves > 0 joined from the raw Wyscout CSV. Recomputes a data-driven
# one-way ICC for all-player vs. outfield-only rows across all 53 features to
# identify which features are confounded by GK structural heterogeneity.
#
# The data-driven ICC (ANOVA decomposition of Z-scaled Y) differs from the
# model-based ICC in 16_icc_summary.csv (posterior variance components from the
# stage-10 Normal fit). Numbers won't match exactly; the comparison of
# icc_all vs. icc_outfield within this script is apples-to-apples and reliable
# for identifying GK-confounded features directionally.
#
# Outputs:
#   outputs/tables/25_gk_icc_comparison.csv
#   outputs/notes/note_goalkeeper_check.md

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/25_gk_icc_analysis.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ---- Load model objects and Y -----------------------------------------------

model_objects <- readRDS(file.path(paths$processed, "model_objects.rds"))
Y <- as.matrix(read_csv(file.path(paths$processed, "Y_scaled.csv"),
                        show_col_types = FALSE))
colnames(Y) <- model_objects$variable_names
N <- nrow(Y)
P <- ncol(Y)
feat_cols <- model_objects$variable_names

meta <- model_objects$metadata |>
  mutate(row_n = row_number())

player_index <- model_objects$player_index   # length N, integer
season_index <- model_objects$season_index

message("Loaded Y: N=", N, " rows, P=", P, " features")

# ---- Build GK flag from raw CSV ---------------------------------------------
# Use total_gk_saves > 0 as the GK indicator. The raw CSV has one row per
# player-match; aggregate to player-season level to get the season total.

message("Reading raw CSV to build GK flag...")
raw <- read_csv(paths$raw,
                col_select  = c(player_id, season_id, total_gk_saves),
                show_col_types = FALSE)

gk_flag <- raw |>
  group_by(player_id, season_id) |>
  summarise(total_gk_saves_season = sum(total_gk_saves, na.rm = TRUE),
            .groups = "drop") |>
  mutate(is_gk = total_gk_saves_season > 0)

# Join onto metadata (meta has one row per model row)
meta_gk <- meta |>
  left_join(gk_flag |> select(player_id, season_id, is_gk),
            by = c("player_id", "season_id")) |>
  mutate(is_gk = replace_na(is_gk, FALSE))

n_gk_rows    <- sum(meta_gk$is_gk)
n_gk_players <- n_distinct(meta_gk$player_id[meta_gk$is_gk])
message(sprintf("GK rows: %d / %d (%.1f%%) across %d distinct GK players",
                n_gk_rows, N, 100 * n_gk_rows / N, n_gk_players))

if (n_gk_rows > 0) {
  message("GK player names (first 20):")
  gk_names <- meta_gk |>
    filter(is_gk) |>
    distinct(player_name, player_id) |>
    arrange(player_name)
  print(head(as.data.frame(gk_names), 20))
}

# ---- One-way random-effects ICC by player -----------------------------------
# ICC = (MS_between - MS_within) / (MS_between - MS_within + n_harmonic * MS_within)
# where n_harmonic = (N - sum(n_i^2)/N) / (I - 1).

compute_icc_player <- function(y, player_id_vec) {
  player_id_vec <- as.character(player_id_vec)
  grand_mean    <- mean(y)
  pm            <- tapply(y, player_id_vec, mean)         # player means
  n_i           <- as.integer(table(player_id_vec))        # obs per player
  N_sub         <- length(y)
  I_sub         <- length(pm)

  if (I_sub <= 1L || N_sub <= I_sub) return(NA_real_)

  # Expand player means to each observation
  pm_per_obs <- pm[player_id_vec]

  SS_between  <- sum(n_i * (pm - grand_mean)^2)
  SS_within   <- sum((y - pm_per_obs)^2)
  MS_between  <- SS_between / (I_sub - 1L)
  MS_within   <- SS_within  / (N_sub - I_sub)
  n_harm      <- (N_sub - sum(n_i^2) / N_sub) / (I_sub - 1L)

  icc <- (MS_between - MS_within) / (MS_between - MS_within + n_harm * MS_within)
  max(0, icc, na.rm = TRUE)
}

# Compute ICC for all rows and outfield-only rows
gk_mask     <- meta_gk$is_gk
outfield_idx <- which(!gk_mask)

message("Computing data-driven ICC for all rows and outfield-only rows...")
icc_all      <- vapply(seq_len(P), function(p)
  compute_icc_player(Y[, p], meta_gk$player_id),
  numeric(1))

icc_outfield <- vapply(seq_len(P), function(p)
  compute_icc_player(Y[outfield_idx, p], meta_gk$player_id[outfield_idx]),
  numeric(1))

result <- tibble(
  feature       = feat_cols,
  icc_all       = round(icc_all, 4),
  icc_outfield  = round(icc_outfield, 4),
  delta         = round(icc_all - icc_outfield, 4)   # positive = GK inflated
)

# Join in the model-based icc_player from stage 16 for reference
icc16 <- read_csv(file.path(paths$tables, "16_icc_summary.csv"),
                  show_col_types = FALSE) |>
  select(feature, icc_player_model = icc_player)

result <- result |>
  left_join(icc16, by = "feature") |>
  arrange(desc(abs(delta)))

write_csv(result, file.path(paths$tables, "25_gk_icc_comparison.csv"))
message("Saved 25_gk_icc_comparison.csv")

# ---- Print summary ----------------------------------------------------------

message("\n--- Features with largest GK-inflation of icc_player (|delta| > 0.02) ---")
confounded <- result |> filter(abs(delta) > 0.02) |> arrange(desc(abs(delta)))
print(as.data.frame(confounded))

message("\n--- Borderline features after GK correction ---")
borderline <- result |>
  filter(icc_player_model >= 0.20, icc_player_model < 0.35) |>
  select(feature, icc_all, icc_outfield, delta, icc_player_model) |>
  arrange(icc_player_model)
print(as.data.frame(borderline))

# ---- Write note file --------------------------------------------------------

note_path <- file.path(paths$notes, "note_goalkeeper_check.md")

gk_confounded_list <- confounded |>
  mutate(line = sprintf("- `%s`: icc_all=%.3f, icc_outfield=%.3f (delta=%.3f)",
                        feature, icc_all, icc_outfield, delta)) |>
  pull(line)

borderline_verdict <- borderline |>
  mutate(verdict = ifelse(icc_outfield < 0.20, "drop after GK fix", "still borderline")) |>
  mutate(line = sprintf("- `%s`: icc_outfield=%.3f → %s",
                        feature, icc_outfield, verdict)) |>
  pull(line)

note_lines <- c(
  "# Goalkeeper Structural Check",
  "",
  "## Finding",
  "",
  sprintf("Goalkeepers were **not** excluded from the model panel. Only GK-specific"),
  sprintf("feature columns (prefix `gk_` or `goal_kicks`) were filtered in stage 01."),
  sprintf("GK flag: `total_gk_saves > 0` joined from the raw Wyscout CSV by player_id/season_id."),
  "",
  sprintf("- **%d GK player-season rows** out of N=%d (%.1f%%)",
          n_gk_rows, N, 100 * n_gk_rows / N),
  sprintf("- **%d distinct GK players** present in `model_objects$metadata`", n_gk_players),
  "",
  "## GK-confounded features (|icc_all - icc_outfield| > 0.02)",
  "",
  paste(gk_confounded_list, collapse = "\n"),
  "",
  "## Impact on borderline features (model-based icc_player 0.20–0.35)",
  "",
  paste(borderline_verdict, collapse = "\n"),
  "",
  "## Recommended fix",
  "",
  "Add a row-level GK filter in `scripts/01_read_clean_feature_engineering.R`",
  "before the `eligible` step, using the `total_gk_saves` column from the",
  "longitudinal panel produced by stage 00 (`data/interim/longitudinal_player_seasons.csv`).",
  "",
  "```r",
  "eligible <- panel |>",
  "  filter(minutes >= pipeline$minimum_minutes_model,",
  "         total_gk_saves == 0 | is.na(total_gk_saves))  # exclude GK rows",
  "```",
  "",
  "This goes in stage 01's rebuild, bundled with Task K's feature drops and transforms.",
  "",
  sprintf("_Written by scripts/25_gk_icc_analysis.R on %s_", Sys.Date())
)

writeLines(note_lines, note_path)
message("Saved note_goalkeeper_check.md")
message("\nTask I2 complete.")
