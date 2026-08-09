# Stage 35 — Minutes-in-variance diagnostic on stage-34a residuals
#
# Re-runs the stage-30 quantile test using minutes-scaled standardised residuals
# from the stage-34a model (fixed φ=0.5). Verifies that the minutes-scaled
# likelihood corrects the low-minute over-representation found at stage 30
# (48.9% of top-200 worst cells from bottom-20th percentile of minutes, p<4.4e-18).
#
# Key difference from stage 30: denominator includes sqrt(m_ref / minutes_n),
# which deflates residuals from low-minute observations.
#
# Column indices for A and B verified against stage-34a CSV header:
#   A[I,P]: cols 157022..230413 (0-indexed) — identical to stage 28
#   B[S,P]: cols 231038..231613 (0-indexed) — identical to stage 28

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/35_minutes_diagnostic_34a.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

# ── Column-range mean extractor (copied from stage 29) ──────────────────────

extract_col_range_means <- function(csv_paths, start_col, end_col) {
  py <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys",
    "paths = sys.argv[1].split(',')",
    "s = int(sys.argv[2]); e = int(sys.argv[3])",
    "width = e - s + 1",
    "acc = [0.0] * width; count = 0",
    "for path in paths:",
    "    with open(path) as f:",
    "        header_done = False",
    "        for line in f:",
    "            if line.startswith('#'): continue",
    "            if not header_done: header_done = True; continue",
    "            vals = line.rstrip().split(',')",
    "            for i in range(width): acc[i] += float(vals[s + i])",
    "            count += 1",
    "for v in acc: print(v / count)"
  ), py)
  on.exit(unlink(py), add = TRUE)
  out <- system2("python3",
                 c(py, paste(csv_paths, collapse = ","),
                   as.character(start_col), as.character(end_col)),
                 stdout = TRUE, stderr = FALSE)
  if (length(out) != (end_col - start_col + 1L))
    stop("extract_col_range_means: expected ", end_col - start_col + 1L,
         " values, got ", length(out))
  as.numeric(out)
}

# ── 1. Load data ─────────────────────────────────────────────────────────────

mo            <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names          # P=48
row_map       <- read_csv(file.path(paths$processed, "model_row_mapping.csv"),
                          show_col_types = FALSE)
Y_raw    <- read_csv(file.path(paths$processed, "Y_scaled.csv"), show_col_types = FALSE)
Y_scaled <- as.matrix(Y_raw[, feature_names])   # N × P

minutes  <- mo$metadata$minutes             # N-vector, 450..3738
m_ref    <- median(minutes)                 # 1585

I <- length(unique(row_map$player_index))   # 1529
S <- length(unique(row_map$season_index))   # 12
P <- length(feature_names)                  # 48
N <- nrow(row_map)                          # 4586

stopifnot(nrow(Y_scaled) == N, ncol(Y_scaled) == P, length(minutes) == N)
message(sprintf("Loaded: N=%d, I=%d, S=%d, P=%d, m_ref=%.0f", N, I, S, P, m_ref))

# ── 2. Stage-34a sigma_e and nu from posterior summary ───────────────────────

summ_path <- file.path(paths$tables,
  "34a_real_lowrank_a_diag_b_t_mv_k3_posterior_summary.csv")
if (!file.exists(summ_path))
  stop("34a posterior summary not found — run stage 34 first.")

summ_34a <- read_csv(summ_path, show_col_types = FALSE)
nu_hat   <- summ_34a$mean[summ_34a$variable == "nu"]
sigma_e  <- summ_34a$mean[grepl("^sigma_e\\[", summ_34a$variable)]

if (length(nu_hat) != 1)  stop("nu not found in 34a posterior summary")
if (length(sigma_e) != P) stop("sigma_e has wrong length — expected ", P, " got ", length(sigma_e))

message(sprintf("Stage-34a: nu=%.4f, sigma_e range [%.4f, %.4f]",
                nu_hat, min(sigma_e), max(sigma_e)))

# Minutes-scaling and residual variance scaling
sqrt_scale_n <- sqrt(m_ref / minutes)              # N-vector, >1 for low-minute
eps_scale    <- sqrt(nu_hat / (nu_hat - 2))        # scalar
tau_99       <- qt(0.995, df = nu_hat) / eps_scale # 99th-pct threshold
message(sprintf("tau_99 at nu=%.4f: %.4f (stage-28 was 3.135)", nu_hat, tau_99))

# ── 3. Extract A and B posterior means from stage-34a CSVs ──────────────────

csv_files <- discover_cmdstan_csv_files("34a_real_lowrank_a_diag_b_t_mv_k3")
if (length(csv_files) == 0L)
  stop("Stage-34a CSV files not found — run stage 34 first.")
message(sprintf("Using %d stage-34a CSV chains for extraction", length(csv_files)))

# Column indices verified against CSV header (same as stage 28 — identical model structure):
#   A[I,P]: cols 157022..230413  (I*P = 73,392 values)
#   B[S,P]: cols 231038..231613  (S*P = 576 values)
message("Extracting A_mean (I×P = 73,392 values) — expect ~30 min ...")
t0 <- proc.time()[["elapsed"]]
A_vec  <- extract_col_range_means(csv_files, 157022L, 230413L)
A_mean <- matrix(A_vec, nrow = I, ncol = P, byrow = FALSE)
colnames(A_mean) <- feature_names
message(sprintf("  A_mean done in %.0f s", proc.time()[["elapsed"]] - t0))

message("Extracting B_mean (S×P = 576 values) ...")
B_vec  <- extract_col_range_means(csv_files, 231038L, 231613L)
B_mean <- matrix(B_vec, nrow = S, ncol = P, byrow = FALSE)
colnames(B_mean) <- feature_names

stopifnot(nrow(A_mean) == I, ncol(A_mean) == P,
          nrow(B_mean) == S, ncol(B_mean) == P)
message(sprintf("A_mean range: [%.3f, %.3f]; B_mean range: [%.3f, %.3f]",
                min(A_mean), max(A_mean), min(B_mean), max(B_mean)))

# ── 4. Minutes-scaled standardised residuals ─────────────────────────────────

A_n     <- A_mean[row_map$player_index, ]    # N × P
B_n     <- B_mean[row_map$season_index, ]    # N × P
eps_raw <- Y_scaled - A_n - B_n             # N × P raw residuals
rm(A_n, B_n)

# z_{n,p} = eps_{n,p} / (sigma_e[p] * sqrt_scale_n * eps_scale)
# Apply: divide by sigma_e[p] (per column), then by sqrt_scale_n * eps_scale (per row)
scale_e_full <- sigma_e * eps_scale                  # P-vector: sigma_e[p] * sqrt(nu/(nu-2))
z_std <- sweep(eps_raw, 2, scale_e_full, "/")        # N × P: divide by sigma_e[p]*eps_scale
z_std <- z_std / sqrt_scale_n                        # N × P: divide by sqrt(m_ref/m_n)
# Equivalent to: z_{n,p} = eps / (sigma_e[p] * sqrt(m_ref/m_n) * sqrt(nu/(nu-2)))

message(sprintf("Residuals (34a): mean=%.4f, sd=%.4f, p99(|z|)=%.4f",
                mean(z_std), sd(z_std), quantile(abs(z_std), 0.99)))
rm(eps_raw)

# ── 5. Top-200 worst cells ───────────────────────────────────────────────────

top_n     <- 200L
abs_z_vec <- as.vector(abs(z_std))              # N*P (column-major)
top_idx   <- order(abs_z_vec, decreasing = TRUE)[seq_len(top_n)]
row_idx   <- ((top_idx - 1L) %% N) + 1L
col_idx   <- ((top_idx - 1L) %/% N) + 1L

worst_tbl <- data.frame(
  rank         = seq_len(top_n),
  player_name  = row_map$player_name[row_idx],
  player_id    = row_map$player_id[row_idx],
  season_id    = row_map$season_id[row_idx],
  feature      = feature_names[col_idx],
  z_std        = z_std[cbind(row_idx, col_idx)],
  minutes      = minutes[row_idx],
  player_index = row_map$player_index[row_idx],
  season_index = row_map$season_index[row_idx]
)

write_csv(worst_tbl, file.path(paths$tables, "35_residual_worst_cells_34a.csv"))
message("Saved: 35_residual_worst_cells_34a.csv")

# ── 6. Quantile share test ───────────────────────────────────────────────────

p20_min <- quantile(minutes, 0.20)   # 889 min (same as stage 30)

obs_worst <- worst_tbl |>
  group_by(player_index, season_index) |>
  summarise(max_abs_z = max(abs(z_std)), minutes = first(minutes), .groups = "drop")

n_total_obs <- nrow(obs_worst)
n_low_min   <- sum(obs_worst$minutes < p20_min, na.rm = TRUE)
share_low   <- n_low_min / n_total_obs

cat(sprintf("\n20th percentile of minutes: %.0f\n", p20_min))
cat(sprintf("Top-200 obs with minutes < p20 (stage-34a): %d / %d = %.1f%%  (expected 20%%)\n",
            n_low_min, n_total_obs, 100 * share_low))
cat(sprintf("Stage-28 baseline was: 48.9%% (88/180 obs)\n"))

# Binomial test: is share_low consistent with 20%?
binom_test <- binom.test(n_low_min, n_total_obs, p = 0.20, alternative = "greater")
cat(sprintf("Binomial test p-value: %.3g\n", binom_test$p.value))

# ── 7. Regression: max|z| ~ log(minutes) ────────────────────────────────────

lm_log  <- lm(max_abs_z ~ log(minutes), data = obs_worst)
smry_log <- summary(lm_log)$coefficients
slope_log <- smry_log["log(minutes)", "Estimate"]
p_log     <- smry_log["log(minutes)", "Pr(>|t|)"]

cat(sprintf("Regression max|z34a| ~ log(min): slope=%.4f (p=%.4f)\n", slope_log, p_log))

# ── 8. Verdict ───────────────────────────────────────────────────────────────

# Full distribution comparison: top-200 vs all observations
all_obs <- data.frame(minutes = minutes)

confirm_slope    <- slope_log < 0 && p_log < 0.05
confirm_quantile <- share_low > 0.30   # >30% still over-represented
corrected        <- share_low < 0.25   # within 5pp of expected 20% → corrected

verdict <- if (corrected) {
  "CORRECTED — minutes-scaled model eliminates low-minute over-representation."
} else if (!confirm_quantile) {
  "LARGELY CORRECTED — share dropped substantially vs stage-28 (48.9%), close to 20% expected."
} else {
  "PARTIAL — share reduced but still over-represented. Consider estimated phi or larger phi."
}

# Compare to stage-28: binomial p at stage 28 was <4.4e-18
# At stage 34a, if corrected: p should be much larger (not significant)

note_lines <- c(
  "=== Stage 35: Minutes-in-variance diagnostic on stage-34a residuals ===",
  sprintf("Date: %s", Sys.Date()),
  "",
  sprintf("Model: Stage 34a (fixed phi=0.5, nu=%.4f)", nu_hat),
  sprintf("Data: N=%d observations, minutes range %.0f-%.0f (median %.0f)",
          N, min(minutes), max(minutes), m_ref),
  "",
  "--- Stage-28 baseline (from stage 30) ---",
  "Low-minute share in top-200 worst: 48.9% (88/180 obs), binomial p < 4.4e-18",
  "",
  "--- Stage-34a (minutes-scaled) ---",
  sprintf("Residual formula: z_{n,p} = (y - A_hat - B_hat) / (sigma_e[p] * sqrt(m_ref/m_n) * sqrt(nu/(nu-2)))"),
  sprintf("tau_99 at nu=%.4f: %.4f", nu_hat, tau_99),
  sprintf("20th pct minutes: %.0f min", p20_min),
  sprintf("Top-200 obs with minutes < p20: %d/%d = %.1f%% (expected 20%%)",
          n_low_min, n_total_obs, 100 * share_low),
  sprintf("Binomial p-value (H1: share > 20%%): %.3g", binom_test$p.value),
  sprintf("Regression max|z| ~ log(min): slope=%.4f, p=%.4f", slope_log, p_log),
  "",
  "--- VERDICT ---",
  verdict,
  "",
  "Figures: 35_minutes_vs_residual_34a.png, 35_minutes_distribution_comparison_34a.png"
)

writeLines(note_lines, file.path(paths$notes, "35_minutes_diagnostic_34a.txt"))
cat(paste(note_lines, collapse = "\n"), "\n")

# ── 9. Figures ───────────────────────────────────────────────────────────────

# Figure 1: max|z34a| vs minutes (per unique obs in top-200)
p_scatter <- ggplot(obs_worst, aes(x = minutes, y = max_abs_z)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "#2166ac") +
  geom_smooth(method = "loess", se = TRUE, colour = "#d73027", linewidth = 0.9) +
  geom_vline(xintercept = p20_min, linetype = "dashed", colour = "grey50") +
  annotate("text", x = p20_min + 30, y = max(obs_worst$max_abs_z) * 0.95,
           hjust = 0, size = 3, colour = "grey40",
           label = sprintf("p20 = %.0f min", p20_min)) +
  labs(
    title    = "Minutes played vs worst standardised residual per player-season (stage 34a)",
    subtitle = sprintf("Minutes-scaled model (phi=0.5, nu=%.3f). Slope log(min): %.3f (p=%.3f)",
                       nu_hat, slope_log, p_log),
    x = "Minutes played (season)",
    y = "max |z_std| across 48 features"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(paths$figures, "35_minutes_vs_residual_34a.png"),
       p_scatter, width = 9, height = 6, dpi = 150)
message("Saved: 35_minutes_vs_residual_34a.png")

# Figure 2: distribution comparison — top-200 vs all
mins_df <- bind_rows(
  data.frame(minutes = all_obs$minutes,   group = sprintf("All observations (N=%d)", N)),
  data.frame(minutes = obs_worst$minutes, group = "In top-200 worst cells (stage-34a)")
)

p_hist <- ggplot(mins_df, aes(x = minutes, fill = group)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = p20_min, linetype = "dashed") +
  scale_fill_manual(
    values = setNames(c("#2166ac", "#d73027"),
                      c(sprintf("All observations (N=%d)", N),
                        "In top-200 worst cells (stage-34a)")),
    name = NULL) +
  labs(
    title    = "Minutes distribution: all vs top-200 worst residuals (stage-34a)",
    subtitle = sprintf("%d / %d (%.0f%%) of worst-obs have minutes < p20 (%.0f min); stage-28: 48.9%%",
                       n_low_min, n_total_obs, 100 * share_low, p20_min),
    x = "Minutes played", y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(paths$figures, "35_minutes_distribution_comparison_34a.png"),
       p_hist, width = 9, height = 5, dpi = 150)
message("Saved: 35_minutes_distribution_comparison_34a.png")

message("Stage 35 complete.")
