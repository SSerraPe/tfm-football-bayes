# Stage 29b — Per-K residual adequacy analysis (Session 14).
#
# K-parameterised variant of scripts/29_outlier_study.R that:
#   * auto-discovers A.[i,p] and B.[j,p] column indices from the CSV header
#     (they differ across K because lambda_a_free size changes)
#   * loads sigma_e[p] and nu_hat from the fit's posterior summary
#   * computes standardised residuals, per-feature tail stats, verdict counts
#
# Selects the fit via BFA_K env var (2, 3, or 4). Maps:
#   K=2 -> 18_real_lowrank_a_diag_b_t
#   K=3 -> 28_real_lowrank_a_diag_b_t_k3
#   K=4 -> 21_real_k4_t
#
# Writes:
#   outputs/tables/29b_k{K}_residual_feature_summary.csv
#   outputs/tables/29b_k{K}_residual_skew_triangulation.csv
#   outputs/tables/29b_kcompare_summary.csv (appended across K)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/29b_kcompare_residuals.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

# ── Config ────────────────────────────────────────────────────────────────────

K <- as.integer(Sys.getenv("BFA_K", unset = "3"))
stopifnot(K %in% c(2L, 3L, 4L))

fit_map <- c(
  "2" = "18_real_lowrank_a_diag_b_t",
  "3" = "28_real_lowrank_a_diag_b_t_k3",
  "4" = "21_real_k4_t"
)
fit_id  <- fit_map[[as.character(K)]]
ps_path <- file.path(paths$tables, paste0(fit_id, "_posterior_summary.csv"))

if (!file.exists(ps_path))
  stop("Posterior summary not found: ", ps_path)

message(sprintf("Stage 29b: K=%d, fit_id=%s", K, fit_id))

# ── Load model objects ────────────────────────────────────────────────────────

mo            <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names
row_map       <- read_csv(file.path(paths$processed, "model_row_mapping.csv"),
                          show_col_types = FALSE)
Y_raw    <- read_csv(file.path(paths$processed, "Y_scaled.csv"), show_col_types = FALSE)
Y_scaled <- as.matrix(Y_raw[, feature_names])

I <- length(unique(row_map$player_index))
S <- length(unique(row_map$season_index))
P <- length(feature_names)
N <- nrow(row_map)

# ── Extract sigma_e and nu from posterior summary ─────────────────────────────

ps <- read_csv(ps_path, show_col_types = FALSE)

sigma_e <- ps |>
  filter(grepl("^sigma_e\\[", variable)) |>
  mutate(idx = as.integer(sub("sigma_e\\[(\\d+)\\]", "\\1", variable))) |>
  arrange(idx) |>
  pull(mean)
stopifnot(length(sigma_e) == P)

nu_hat <- ps |> filter(variable == "nu") |> pull(mean)
stopifnot(length(nu_hat) == 1L)

message(sprintf("Loaded: N=%d, I=%d, S=%d, P=%d, K=%d, nu=%.4f", N, I, S, P, K, nu_hat))

# ── Discover A/B column indices from CSV header ───────────────────────────────

csv_files <- discover_cmdstan_csv_files(fit_id)
if (length(csv_files) == 0L) stop("No CSV chains for fit_id=", fit_id)
message(sprintf("Using %d CSV chains for extraction", length(csv_files)))

# Read the first CSV's header (first non-comment line)
first_csv <- csv_files[1]
con <- file(first_csv, "r")
on.exit(close(con), add = TRUE)
header <- NULL
repeat {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (length(line) == 0L) stop("Reached EOF before header")
  if (!startsWith(line, "#")) { header <- line; break }
}
close(con); on.exit()
cols <- strsplit(header, ",", fixed = TRUE)[[1]]

# 0-indexed column positions
find_index <- function(pattern) which(cols == pattern) - 1L
a_first <- find_index("A.1.1")
a_last  <- find_index(sprintf("A.%d.%d", I, P))
b_first <- find_index("B.1.1")
b_last  <- find_index(sprintf("B.%d.%d", S, P))
if (length(a_first) == 0L || length(a_last) == 0L ||
    length(b_first) == 0L || length(b_last) == 0L)
  stop("Could not locate A/B parameter columns in the CSV header of ", first_csv)

# Verify contiguous
if (a_last - a_first + 1L != I * P)
  stop(sprintf("A column span mismatch: %d cols expected, %d found",
               I * P, a_last - a_first + 1L))
if (b_last - b_first + 1L != S * P)
  stop(sprintf("B column span mismatch: %d cols expected, %d found",
               S * P, b_last - b_first + 1L))

message(sprintf("A block cols [%d..%d] (%d values), B block cols [%d..%d] (%d values)",
                a_first, a_last, I * P, b_first, b_last, S * P))

# ── Streaming Python column-range extractor (identical to stage 29) ───────────

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
                 c(py,
                   paste(csv_paths, collapse = ","),
                   as.character(start_col),
                   as.character(end_col)),
                 stdout = TRUE, stderr = FALSE)
  if (length(out) != (end_col - start_col + 1L))
    stop("extract_col_range_means: expected ", end_col - start_col + 1L,
         " values, got ", length(out))
  as.numeric(out)
}

# ── Extract A_mean[I,P] and B_mean[S,P] ───────────────────────────────────────

message("Extracting A_mean (player effects, I×P = ", I * P, " values) ...")
t0 <- proc.time()[["elapsed"]]
A_vec  <- extract_col_range_means(csv_files, a_first, a_last)
A_mean <- matrix(A_vec, nrow = I, ncol = P, byrow = FALSE)
colnames(A_mean) <- feature_names
message(sprintf("  done in %.0f s", proc.time()[["elapsed"]] - t0))

message("Extracting B_mean (season effects, S×P = ", S * P, " values) ...")
B_vec  <- extract_col_range_means(csv_files, b_first, b_last)
B_mean <- matrix(B_vec, nrow = S, ncol = P, byrow = FALSE)
colnames(B_mean) <- feature_names

# Sanity checks
if (abs(mean(A_mean)) > 0.05)
  warning("mean(A_mean) = ", round(mean(A_mean), 4), " — player effects not centred near 0")
if (max(abs(B_mean)) > 5)
  warning("max|B_mean| = ", round(max(abs(B_mean)), 2), " — season effects larger than expected")
message(sprintf("A_mean range: [%.3f, %.3f]; B_mean range: [%.3f, %.3f]",
                min(A_mean), max(A_mean), min(B_mean), max(B_mean)))

# ── Standardised residuals ────────────────────────────────────────────────────

scale_e <- sigma_e * sqrt(nu_hat / (nu_hat - 2))
tau_99  <- qt(0.995, df = nu_hat) / sqrt(nu_hat / (nu_hat - 2))
message(sprintf("tau_99 at nu=%.4f: %.4f", nu_hat, tau_99))

A_n     <- A_mean[row_map$player_index, ]
B_n     <- B_mean[row_map$season_index, ]
eps_raw <- Y_scaled - A_n - B_n
z_std   <- sweep(eps_raw, 2, scale_e, "/")
rm(A_n, B_n)

message(sprintf("Residuals: mean=%.4f, sd=%.4f, p99(|z|)=%.4f",
                mean(z_std), sd(z_std), quantile(abs(z_std), 0.99)))

# ── Per-feature summary ───────────────────────────────────────────────────────

feat_smry <- tibble(feature = feature_names) |>
  mutate(
    p50_abs       = vapply(seq_len(P), function(p) median(abs(z_std[, p])),           numeric(1)),
    p90_abs       = vapply(seq_len(P), function(p) quantile(abs(z_std[, p]), 0.90), numeric(1)),
    p95_abs       = vapply(seq_len(P), function(p) quantile(abs(z_std[, p]), 0.95), numeric(1)),
    p99_abs       = vapply(seq_len(P), function(p) quantile(abs(z_std[, p]), 0.99), numeric(1)),
    max_abs       = vapply(seq_len(P), function(p) max(abs(z_std[, p])),               numeric(1)),
    frac_gt_tau99 = vapply(seq_len(P), function(p) mean(abs(z_std[, p]) > tau_99),  numeric(1))
  ) |>
  arrange(desc(p99_abs))

write_csv(feat_smry, file.path(paths$tables,
          sprintf("29b_k%d_residual_feature_summary.csv", K)))
message(sprintf("Saved: 29b_k%d_residual_feature_summary.csv", K))

# ── Verdict classification (mirrors stage 29's rule) ──────────────────────────

skew_tbl <- read_csv(file.path(paths$tables, "23_feature_skewness.csv"),
                     show_col_types = FALSE) |>
  filter(!is.na(feature))
icc_tbl  <- read_csv(file.path(paths$tables, "16_icc_summary.csv"),
                     show_col_types = FALSE) |>
  filter(!is.na(feature))

# Season concentration: for worst 200 cells, how many distinct seasons?
top_n     <- 200L
abs_z_vec <- as.vector(abs(z_std))
top_idx   <- order(abs_z_vec, decreasing = TRUE)[seq_len(top_n)]
row_idx   <- ((top_idx - 1L) %% N) + 1L
col_idx   <- ((top_idx - 1L) %/% N) + 1L

worst_tbl <- tibble(
  feature   = feature_names[col_idx],
  season_id = row_map$season_id[row_idx]
)
season_conc <- worst_tbl |>
  group_by(feature) |>
  slice_head(n = 20) |>
  summarise(n_distinct_seasons = n_distinct(season_id), .groups = "drop")

tri <- feat_smry |>
  left_join(skew_tbl |> select(feature, skewness, recommended_transform), by = "feature") |>
  left_join(icc_tbl  |> select(feature, icc_player, icc_season),          by = "feature") |>
  left_join(season_conc, by = "feature") |>
  mutate(
    tau_99_val = tau_99,
    verdict = case_when(
      p99_abs <= tau_99                                          ~ "unremarkable_p99",
      icc_player >= 0.35 & p99_abs > tau_99                       ~ "contamination-driven",
      icc_season >= 0.20 & !is.na(n_distinct_seasons) &
        n_distinct_seasons <= 2                                    ~ "season-driven",
      !is.na(recommended_transform) &
        recommended_transform != "none" & p99_abs > tau_99        ~ "skew-residual",
      TRUE                                                         ~ "mixed"
    )
  )
write_csv(tri, file.path(paths$tables,
          sprintf("29b_k%d_residual_skew_triangulation.csv", K)))
message(sprintf("Saved: 29b_k%d_residual_skew_triangulation.csv", K))

# ── Cross-K comparison summary row ────────────────────────────────────────────

summary_row <- tibble(
  K                             = K,
  fit_id                        = fit_id,
  nu_hat                        = nu_hat,
  tau_99                        = tau_99,
  mean_frac_gt_tau99            = mean(feat_smry$frac_gt_tau99),
  median_frac_gt_tau99          = median(feat_smry$frac_gt_tau99),
  max_frac_gt_tau99             = max(feat_smry$frac_gt_tau99),
  expected_frac_under_t         = 0.01,
  n_unremarkable                = sum(tri$verdict == "unremarkable_p99"),
  n_contamination_driven        = sum(tri$verdict == "contamination-driven"),
  n_season_driven               = sum(tri$verdict == "season-driven"),
  n_skew_residual               = sum(tri$verdict == "skew-residual"),
  n_mixed                       = sum(tri$verdict == "mixed"),
  worst_cell_zstd               = max(abs(z_std))
)

comp_path <- file.path(paths$tables, "29b_kcompare_summary.csv")
existing <- if (file.exists(comp_path)) {
  read_csv(comp_path, show_col_types = FALSE) |> filter(K != !!K)
} else NULL
bind_rows(existing, summary_row) |> arrange(K) |> write_csv(comp_path)
message("Appended into 29b_kcompare_summary.csv")

cat("\n=== Cross-K summary row (K=", K, ") ===\n", sep = "")
print(summary_row)

message(sprintf("\nStage 29b (K=%d) complete.", K))
