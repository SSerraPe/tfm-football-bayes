# Stage 22: Residual extraction and feature ICC screen.
#
# TASK A — Residual diagnostics from the stage-18 t-fit.
#   Computes ε_n = y_n − (â_i + b̂_j) using posterior means of player and season
#   effects. Player/season effect means come from stages 14 and 15 (stage-10 Normal
#   fit); these are a very close approximation to the stage-18 t-model means because
#   the mean structure (a_i, b_j) differs only slightly between Normal and t.
#   Standardisation uses sigma_e and nu pulled directly from the stage-18 fit so
#   the t-distribution's Var(ε) = sigma_e² × ν/(ν−2) is applied correctly.
#   Outputs the 200 worst-fitting (player-season, feature) cells and a per-feature
#   tail summary to identify which features carry the heaviest residuals.
#
# TASK B — Feature ICC screen.
#   Reads the existing 16_icc_summary.csv and classifies rate_* and flagged per90_*
#   features into keep / drop_candidate / borderline buckets.
#   NOTE: this is a SCREEN, not a decision. Removing any feature requires a full
#   rebuild from stage 02 onward and must be approved before implementation.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/22_residual_analysis.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ============================================================
# TASK A — Residual extraction
# ============================================================

message("=== Task A: Residual extraction from stage-18 t-fit ===")

# ---- Load Y and index vectors from model_objects --------------------------------

model_objects <- readRDS(file.path(paths$processed, "model_objects.rds"))
Y             <- as.matrix(read_csv(file.path(paths$processed, "Y_scaled.csv"),
                                    show_col_types = FALSE))
colnames(Y)   <- model_objects$variable_names

player_idx  <- model_objects$player_index    # length N; maps row → player (1..I)
season_idx  <- model_objects$season_index    # length N; maps row → season (1..S)
N  <- nrow(Y)
P  <- ncol(Y)

# Build a lookup: model row → player_name + season_id
row_meta <- model_objects$row_mapping |>
  select(model_row = 1, player_index, season_index) |>
  left_join(model_objects$metadata |>
              select(player_id, player_name) |>
              distinct(),
            by = character())  # join will be done below via player_index

# Simpler: use metadata directly (same row order as Y)
meta <- model_objects$metadata |>
  mutate(row_n = row_number()) |>
  select(row_n, player_id, player_name, season_id)

# ---- Load player effect means (stage 14 = Normal fit; close approx to t-fit) ---

a_wide <- read_csv(file.path(paths$tables, "14_player_effect_means.csv"),
                   show_col_types = FALSE)
# Columns: player_id, <feature1>, ..., <featureP>, player_index
feat_cols <- model_objects$variable_names
A_mean_mat <- as.matrix(a_wide[, feat_cols])   # rows = players, cols = features
rownames(A_mean_mat) <- a_wide$player_index     # 1-based player index

# ---- Load season effect means (stage 15 = Normal fit) --------------------------

b_long <- read_csv(file.path(paths$tables, "15_season_effect_means.csv"),
                   show_col_types = FALSE)
# Columns: season, feature, b_mean, ...  where 'season' is the Wyscout season_id
season_map <- model_objects$season_map   # tibble: season_index, season_id
b_long <- b_long |>
  left_join(season_map, by = c("season" = "season_id")) |>
  filter(!is.na(season_index))

B_mean_mat <- b_long |>
  pivot_wider(id_cols = season_index, names_from = feature, values_from = b_mean) |>
  arrange(season_index)
B_mean_mat <- as.matrix(B_mean_mat[, feat_cols])   # rows = seasons, cols = features
rownames(B_mean_mat) <- sort(unique(b_long$season_index))

# ---- Load sigma_e and nu from pre-extracted cache (NOT from ICC table) ---------
# The full stage-18 CmdStan CSV files are ~12GB and loading via as_cmdstan_fit()
# parses all 276,606 columns × 4,000 draws, taking 45+ minutes. Instead we use
# a pre-extracted cache written by a Python awk-style pass over the CSV columns
# (nu = col 8, sigma_e.1..53 = cols 91606–91658 in the CmdStan CSV).
# The cache holds posterior means across all 4 chains × 1000 draws.

cache_path <- file.path(paths$tables, "18_sigma_e_nu_cache.csv")
if (!file.exists(cache_path)) {
  stop("18_sigma_e_nu_cache.csv not found. Run the Python cache-extraction step first.")
}

cache <- read_csv(cache_path, show_col_types = FALSE)
# Ensure order matches feat_cols (should already match, but be explicit)
cache <- cache[match(feat_cols, cache$feature), ]

sigma_e_mean <- cache$sigma_e_mean   # length P; posterior mean from stage-18 t-fit
nu_mean      <- cache$nu_mean[1]     # scalar ≈ 2.95

# t-distribution residual variance: Var(ε) = sigma_e² × ν/(ν−2)
var_e   <- sigma_e_mean^2 * nu_mean / (nu_mean - 2)
sd_e    <- sqrt(var_e)

message(sprintf("nu (stage-18 posterior mean): %.4f", nu_mean))
message(sprintf("sigma_e range: [%.3f, %.3f]; sd_e (t-scaled) range: [%.3f, %.3f]",
                min(sigma_e_mean), max(sigma_e_mean), min(sd_e), max(sd_e)))

# ---- Compute raw and standardised residuals ------------------------------------

message("Computing residuals for N=", N, " observations × P=", P, " features...")

# eps_raw[n, p] = Y[n, p] - A_mean[player_idx[n], p] - B_mean[season_idx[n], p]
A_row <- A_mean_mat[player_idx, ]   # N × P; each row is the player's mean effect
B_row <- B_mean_mat[season_idx, ]   # N × P; each row is the season's mean effect

eps_raw <- Y - A_row - B_row        # N × P raw residuals

# Standardise: eps_std[n, p] = eps_raw[n, p] / sd_e[p]
eps_std <- sweep(eps_raw, 2, sd_e, "/")   # N × P standardised residuals

# ---- Rank worst-fitting cells ---------------------------------------------------

message("Identifying worst-fitting cells...")

abs_eps <- abs(eps_std)

# Long format: one row per (n, p) observation
long <- as.data.frame(which(!is.na(abs_eps), arr.ind = TRUE))
colnames(long) <- c("row_n", "feat_idx")
long$feature   <- feat_cols[long$feat_idx]
long$eps_raw   <- eps_raw[cbind(long$row_n, long$feat_idx)]
long$eps_std   <- eps_std[cbind(long$row_n, long$feat_idx)]
long$abs_std   <- abs(long$eps_std)

long <- long |>
  left_join(meta, by = "row_n") |>
  left_join(model_objects$season_map, by = c("season_id" = "season_id")) |>
  arrange(desc(abs_std))

worst200 <- long |>
  select(player_name, player_id, season_id, season_index, feature,
         eps_std, eps_raw, abs_std) |>
  head(200)

write_csv(worst200, file.path(paths$tables, "22_residual_worst_cells.csv"))
message("Saved 22_residual_worst_cells.csv (top 200 worst-fitting cells)")

# ---- Per-feature tail summary ---------------------------------------------------

feat_summary <- long |>
  group_by(feature) |>
  summarise(
    median_abs_std = median(abs_std),
    p90            = quantile(abs_std, 0.90),
    p95            = quantile(abs_std, 0.95),
    p99            = quantile(abs_std, 0.99),
    max_abs_std    = max(abs_std),
    n_gt3sd        = sum(abs_std > 3),
    n_gt5sd        = sum(abs_std > 5),
    .groups = "drop"
  ) |>
  arrange(desc(p99))

write_csv(feat_summary, file.path(paths$tables, "22_residual_feature_summary.csv"))
message("Saved 22_residual_feature_summary.csv (per-feature tail summary)")

# Print summary
message("\nTop 10 features by |eps_std| p99 (heaviest-tail features):")
print(as.data.frame(head(feat_summary, 10)))

message("\nTop 20 worst-fitting cells:")
print(as.data.frame(head(worst200 |> select(player_name, season_id, feature, eps_std), 20)))

# ============================================================
# TASK B — Feature ICC screen
# ============================================================

message("\n=== Task B: Feature ICC screen ===")

icc <- read_csv(file.path(paths$tables, "16_icc_summary.csv"),
                show_col_types = FALSE)

# Classification rules:
#   icc_player < 0.20  → drop_candidate (nearly all residual noise)
#   0.20 – 0.35        → borderline
#   > 0.35             → keep
# High icc_season (> 0.30) is flagged separately as "season_driven" even if
# icc_player is moderate (signal is real but not stable player identity).

classify_feature <- function(icc_player, icc_season, feature) {
  verdict <- dplyr::case_when(
    icc_player < 0.20              ~ "drop_candidate",
    icc_player < 0.35              ~ "borderline",
    TRUE                           ~ "keep"
  )
  season_flag <- ifelse(icc_season > 0.30, "season_driven", "")
  paste(verdict, season_flag, sep = ifelse(season_flag != "", ";", ""))
}

icc_screen <- icc |>
  mutate(
    feature_type = ifelse(grepl("^rate_", feature), "rate", "per90"),
    verdict = classify_feature(icc_player, icc_season, feature)
  ) |>
  select(feature, feature_type, group, icc_player, icc_season, icc_residual, verdict) |>
  arrange(icc_player)

write_csv(icc_screen, file.path(paths$tables, "22_icc_screen_summary.csv"))
message("Saved 22_icc_screen_summary.csv")

message("\n--- rate_* features (sorted by icc_player) ---")
print(as.data.frame(icc_screen |> filter(feature_type == "rate") |>
                      select(feature, icc_player, icc_residual, verdict)))

message("\n--- per90_* features with icc_player < 0.50 ---")
print(as.data.frame(icc_screen |>
                      filter(feature_type == "per90", icc_player < 0.50) |>
                      select(feature, icc_player, icc_season, verdict)))

message("\nACTION NOTE: Removing any feature requires stage 02+ rebuild. This is a",
        " screen only — present to professor before dropping anything.")

message("\nStage 22 complete.")
