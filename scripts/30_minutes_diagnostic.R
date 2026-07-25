# Stage 30 — Minutes-in-variance diagnostic (Task 1a, Session 11)
#
# Hypothesis (Víctor): per90 features are ratio estimators; a player-season with few
# minutes has higher measurement variance ∝ 1/minutes. The current likelihood uses
# constant σ_{e,p}, so low-minute rows may look heavy-tailed — driving ν low and
# dominating the outlier list.
#
# This script gates Task 1b (minutes-scaled Stan model). No refit needed.
# Verdict written to outputs/notes/30_minutes_diagnostic.txt.

source("src/bootstrap.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# ── 1. Load data ──────────────────────────────────────────────────────────────

mo      <- readRDS(file.path(paths$processed, "model_objects.rds"))
meta    <- mo$metadata                           # N=4586 rows, indexed by observation n
minutes <- meta$minutes                          # N-vector, 450–3738

cat(sprintf("N observations: %d\n", nrow(meta)))
cat(sprintf("Minutes range: %.0f – %.0f  (median %.0f)\n",
            min(minutes), max(minutes), median(minutes)))
cat(sprintf("20th pct minutes: %.0f\n", quantile(minutes, 0.20)))

# ── 2. Per-observation worst |z| from stage-29 worst-cells ────────────────────
# 29_residual_worst_cells.csv has top-200 cells; group by (player_index, season_index)
# to get the worst |z| seen for each observation.

wc_path <- file.path(paths$tables, "29_residual_worst_cells.csv")
if (!file.exists(wc_path)) stop("29_residual_worst_cells.csv not found — run stage 29 first.")

wc <- read_csv(wc_path, show_col_types = FALSE)

# Top |z| per observation (player_index × season_index identifies a unique row_n)
obs_worst <- wc |>
  group_by(player_index, season_index) |>
  summarise(max_abs_z = max(abs(z_std)), .groups = "drop")

cat(sprintf("Distinct observations in worst-200: %d\n", nrow(obs_worst)))

# Join minutes: meta is ordered by observation n (1:N).
# row_map maps n → player_index, season_index
row_map_path <- file.path(paths$processed, "model_row_mapping.csv")
row_map <- read_csv(row_map_path, show_col_types = FALSE)

meta_full <- row_map |>
  mutate(minutes = minutes[seq_len(n())]) |>  # meta is same order as row_map
  select(player_index, season_index, minutes)

obs_worst <- obs_worst |>
  left_join(meta_full, by = c("player_index", "season_index"))

cat(sprintf("Rows with minutes after join: %d  (NA: %d)\n",
            sum(!is.na(obs_worst$minutes)), sum(is.na(obs_worst$minutes))))

# Also build the full observation-level minutes table (all N rows)
all_obs <- meta_full

# ── 3. Regressions ────────────────────────────────────────────────────────────

dat <- obs_worst |> filter(!is.na(minutes), !is.na(max_abs_z))

lm_log  <- lm(max_abs_z ~ log(minutes),    data = dat)
lm_inv  <- lm(max_abs_z ~ I(1 / minutes),  data = dat)

smry_log <- summary(lm_log)$coefficients
smry_inv <- summary(lm_inv)$coefficients

cat("\n=== Regression: max|z| ~ log(minutes) ===\n")
print(round(smry_log, 4))

cat("\n=== Regression: max|z| ~ 1/minutes ===\n")
print(round(smry_inv, 4))

slope_log <- smry_log["log(minutes)", "Estimate"]
p_log     <- smry_log["log(minutes)", "Pr(>|t|)"]
slope_inv <- smry_inv["I(1/minutes)", "Estimate"]
p_inv     <- smry_inv["I(1/minutes)", "Pr(>|t|)"]

# ── 4. Quantile share: how many of the top-200 worst cells are low-minute? ───

p20_minutes <- quantile(all_obs$minutes, 0.20)

# All N observations in top-200 (one per unique player_index×season_index)
n_total_obs <- nrow(obs_worst)
n_low_min   <- sum(obs_worst$minutes < p20_minutes, na.rm = TRUE)
share_low   <- n_low_min / n_total_obs
expected    <- 0.20   # under uniform distribution

cat(sprintf("\n20th percentile of minutes: %.0f\n", p20_minutes))
cat(sprintf("Obs in worst-200 with minutes < p20: %d / %d = %.1f%%  (expected under uniform: 20%%)\n",
            n_low_min, n_total_obs, 100 * share_low))

# ── 5. Figures ────────────────────────────────────────────────────────────────

# Figure 1: scatter max|z| vs minutes with LOESS
p_scatter <- ggplot(dat, aes(x = minutes, y = max_abs_z)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#2166ac") +
  geom_smooth(method = "loess", se = TRUE, colour = "#d73027", linewidth = 0.9) +
  geom_vline(xintercept = p20_minutes, linetype = "dashed", colour = "grey50") +
  annotate("text", x = p20_minutes + 30, y = max(dat$max_abs_z, na.rm = TRUE) * 0.95,
           hjust = 0, size = 3, colour = "grey40",
           label = sprintf("p20 = %.0f min", p20_minutes)) +
  labs(
    title    = "Minutes played vs worst standardised residual per player-season",
    subtitle = sprintf("Stage 28 (K=3, t-model, ν=4.902). Slope log(min): %.3f (p=%.3f)",
                       slope_log, p_log),
    x = "Minutes played (season)",
    y = "max |z_std| across 48 features"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(paths$figures, "30_minutes_vs_residual.png"),
       p_scatter, width = 9, height = 6, dpi = 150)
message("Saved: 30_minutes_vs_residual.png")

# Figure 2: minutes distribution — top-200 obs vs all observations
minutes_df <- bind_rows(
  data.frame(minutes = all_obs$minutes,     group = "All observations (N=4586)"),
  data.frame(minutes = obs_worst$minutes,   group = "In top-200 worst cells")
)

p_hist <- ggplot(minutes_df, aes(x = minutes, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = p20_minutes, linetype = "dashed") +
  scale_fill_manual(values = c("All observations (N=4586)" = "#2166ac",
                                "In top-200 worst cells"    = "#d73027"),
                    name = NULL) +
  labs(
    title    = "Minutes distribution: all observations vs top-200 worst residuals",
    subtitle = sprintf("%d / %d (%.0f%%) of worst-obs have minutes < p20 (%.0f min); expected 20%%",
                       n_low_min, n_total_obs, 100 * share_low, p20_minutes),
    x = "Minutes played", y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(paths$figures, "30_minutes_distribution_comparison.png"),
       p_hist, width = 9, height = 5, dpi = 150)
message("Saved: 30_minutes_distribution_comparison.png")

# ── 6. Verdict ────────────────────────────────────────────────────────────────

confirm_slope    <- slope_log < 0 && p_log < 0.05
confirm_quantile <- share_low > 0.30   # >30% in bottom 20% → substantially over-represented

verdict <- if (confirm_slope && confirm_quantile) {
  "CONFIRM — negative slope + over-represented low-minute obs. Proceed with Task 1b."
} else if (confirm_slope || confirm_quantile) {
  "PARTIAL — one signal positive. Investigate further before committing to Task 1b."
} else {
  "REJECT — no clear minutes gradient. Task 1b not warranted."
}

note_lines <- c(
  "=== Stage 30: Minutes-in-variance diagnostic ===",
  sprintf("Date: %s", Sys.Date()),
  "",
  sprintf("Data: N=%d observations, minutes range %.0f–%.0f (median %.0f)",
          nrow(meta), min(minutes), max(minutes), median(minutes)),
  "",
  "--- Regressions ---",
  sprintf("max|z| ~ log(minutes):   slope=%.4f, SE=%.4f, p=%.4f",
          slope_log, smry_log["log(minutes)", "Std. Error"], p_log),
  sprintf("max|z| ~ 1/minutes:      slope=%.4f, SE=%.4f, p=%.4f",
          slope_inv, smry_inv["I(1/minutes)", "Std. Error"], p_inv),
  "",
  "--- Quantile share ---",
  sprintf("20th pct minutes: %.0f min", p20_minutes),
  sprintf("Top-200 obs with minutes < p20: %d/%d = %.1f%% (expected 20%%)",
          n_low_min, n_total_obs, 100 * share_low),
  "",
  "--- VERDICT ---",
  verdict,
  "",
  "Figures: 30_minutes_vs_residual.png, 30_minutes_distribution_comparison.png"
)

writeLines(note_lines, file.path(paths$notes, "30_minutes_diagnostic.txt"))
cat(paste(note_lines, collapse = "\n"), "\n")
message("Stage 30 complete.")
