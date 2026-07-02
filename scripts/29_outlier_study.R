# Stage 29 — Outlier Study on the Final Model (Stage 28, K=3, P=48, N=4586, nu≈4.902)
#
# Extracts posterior-mean player and season effects from stage 28 CSVs via a streaming
# Python column-range reader (avoids loading 9.6 GB into R). Computes standardised
# residuals z[n,p] ~ t(nu)/sqrt(nu/(nu-2)) with unit variance. Produces:
#   29_residual_feature_summary.csv   — per-feature tail statistics
#   29_residual_worst_cells.csv       — top 200 (player-season, feature) cells by |z|
#   29_residual_skew_triangulation.csv — feature classification
#   29_player_factor_scores.csv       — PC scores for all 1529 players
#   29_top_players_by_pc.csv          — top / bottom 20 per PC
#   29_season_factor_trends.csv       — season-level PC projection of B_mean
#   figures: 29_residual_p99_barchart.png, 29_factor_score_scatter_pc1_pc2.png,
#            29_top_players_table_pc{1,2,3}.png, 29_season_factor_trends.png

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/29_outlier_study.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "tidyr"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(forcats)
})

# ── Streaming column-range mean extractor ────────────────────────────────────
# Reads CmdStan CSV files (save_warmup=FALSE — all data rows are sampling draws).
# Skips comment lines and header; computes column means for 0-indexed [start_col, end_col].
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

# ── Step 1 — Data loading ─────────────────────────────────────────────────────

mo            <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names          # length 48 (in config.R row order)
row_map       <- read_csv(file.path(paths$processed, "model_row_mapping.csv"),
                          show_col_types = FALSE)
# row_map: model_row, player_index, season_index, player_id, player_name, season_id

Y_raw <- read_csv(file.path(paths$processed, "Y_scaled.csv"), show_col_types = FALSE)
Y_scaled <- as.matrix(Y_raw[, feature_names])   # N=4586 × P=48, rows = model_row order

dc      <- read_csv(file.path(paths$tables, "28_diagnostics_cache.csv"),
                    show_col_types = FALSE)
sigma_e <- dc |>
  filter(grepl("^sigma_e\\.", param)) |>
  mutate(idx = as.integer(sub("sigma_e\\.(\\d+)", "\\1", param))) |>
  arrange(idx) |>
  pull(mean)                                   # length P=48
nu_hat  <- dc |> filter(param == "nu") |> pull(mean)  # 4.902

I <- length(unique(row_map$player_index))    # 1529
S <- length(unique(row_map$season_index))    # 12
P <- length(feature_names)                   # 48
N <- nrow(row_map)                           # 4586
K <- 3L

stopifnot(nrow(Y_scaled) == N, ncol(Y_scaled) == P, length(sigma_e) == P)
message(sprintf("Loaded: N=%d, I=%d, S=%d, P=%d, K=%d, nu=%.4f", N, I, S, P, K, nu_hat))

# ── Step 2 — Extract A_mean[I,P] and B_mean[S,P] from stage 28 CSVs ──────────

csv_files <- discover_cmdstan_csv_files("28_real_lowrank_a_diag_b_t_k3")
message(sprintf("Using %d stage-28 CSV chains for extraction", length(csv_files)))

# A.i.p: cols 157022-230413 (0-indexed, verified via header probe)
# Stan stores matrix[I,P] column-major: block k = A[1..I, k] (k = feature index)
message("Extracting A_mean (player effects, I×P = 73,392 values) ...")
t0 <- proc.time()[["elapsed"]]
A_vec  <- extract_col_range_means(csv_files, 157022L, 230413L)
A_mean <- matrix(A_vec, nrow = I, ncol = P, byrow = FALSE)  # column-major → correct
colnames(A_mean) <- feature_names
message(sprintf("  done in %.0f s", proc.time()[["elapsed"]] - t0))

# B.j.p: cols 231038-231613 (0-indexed)
message("Extracting B_mean (season effects, S×P = 576 values) ...")
B_vec  <- extract_col_range_means(csv_files, 231038L, 231613L)
B_mean <- matrix(B_vec, nrow = S, ncol = P, byrow = FALSE)
colnames(B_mean) <- feature_names

# Sanity checks
stopifnot(nrow(A_mean) == I, ncol(A_mean) == P,
          nrow(B_mean) == S, ncol(B_mean) == P)
if (abs(mean(A_mean)) > 0.05)
  warning("mean(A_mean) = ", round(mean(A_mean), 4), " — player effects not centred near 0")
if (max(abs(B_mean)) > 5)
  warning("max|B_mean| = ", round(max(abs(B_mean)), 2), " — season effects larger than expected")
message(sprintf("A_mean range: [%.3f, %.3f]; B_mean range: [%.3f, %.3f]",
                min(A_mean), max(A_mean), min(B_mean), max(B_mean)))

# ── Step 3 — Standardised residuals ──────────────────────────────────────────

scale_e <- sigma_e * sqrt(nu_hat / (nu_hat - 2))
tau_99  <- qt(0.995, df = nu_hat) / sqrt(nu_hat / (nu_hat - 2))
message(sprintf("tau_99 at nu=%.4f: %.4f", nu_hat, tau_99))  # expect ~3.135

# Vectorised: broadcast A_mean and B_mean to N rows using player/season indices
A_n     <- A_mean[row_map$player_index, ]   # N × P
B_n     <- B_mean[row_map$season_index, ]   # N × P
eps_raw <- Y_scaled - A_n - B_n             # N × P raw residuals
z_std   <- sweep(eps_raw, 2, scale_e, "/")  # N × P standardised residuals
rm(A_n, B_n)

message(sprintf("Residuals: mean=%.4f, sd=%.4f, p99(|z|)=%.4f",
                mean(z_std), sd(z_std), quantile(abs(z_std), 0.99)))

# ── Step 4 — Feature-level summary ───────────────────────────────────────────

groups <- feature_group_lookup()

feat_smry <- tibble(feature = feature_names) |>
  left_join(groups, by = "feature") |>
  mutate(group = coalesce(group, "Other")) |>
  mutate(
    p50_abs       = vapply(seq_len(P), function(p) median(abs(z_std[, p])),     numeric(1)),
    p90_abs       = vapply(seq_len(P), function(p) quantile(abs(z_std[, p]), 0.90), numeric(1)),
    p95_abs       = vapply(seq_len(P), function(p) quantile(abs(z_std[, p]), 0.95), numeric(1)),
    p99_abs       = vapply(seq_len(P), function(p) quantile(abs(z_std[, p]), 0.99), numeric(1)),
    max_abs       = vapply(seq_len(P), function(p) max(abs(z_std[, p])),         numeric(1)),
    frac_gt_tau99 = vapply(seq_len(P), function(p) mean(abs(z_std[, p]) > tau_99), numeric(1))
  ) |>
  arrange(desc(p99_abs))

write_csv(feat_smry, file.path(paths$tables, "29_residual_feature_summary.csv"))
message("Saved: 29_residual_feature_summary.csv")

# ── Step 5 — Worst-cell ranking ───────────────────────────────────────────────

top_n <- 200L
abs_z_vec <- as.vector(abs(z_std))               # N*P vector (column-major)
top_idx   <- order(abs_z_vec, decreasing = TRUE)[seq_len(top_n)]
row_idx   <- ((top_idx - 1L) %% N) + 1L          # 1..N (col-major: row varies fastest)
col_idx   <- ((top_idx - 1L) %/% N) + 1L         # 1..P

worst_tbl <- tibble(
  rank         = seq_len(top_n),
  player_name  = row_map$player_name[row_idx],
  player_id    = row_map$player_id[row_idx],
  season_id    = row_map$season_id[row_idx],
  feature      = feature_names[col_idx],
  z_std        = z_std[cbind(row_idx, col_idx)],
  eps_raw      = eps_raw[cbind(row_idx, col_idx)],
  player_index = row_map$player_index[row_idx],
  season_index = row_map$season_index[row_idx]
)

write_csv(worst_tbl, file.path(paths$tables, "29_residual_worst_cells.csv"))
message("Saved: 29_residual_worst_cells.csv")

# ── Step 6 — Feature classification ──────────────────────────────────────────

skew_tbl <- read_csv(file.path(paths$tables, "23_feature_skewness.csv"),
                     show_col_types = FALSE) |>
  filter(!is.na(feature))
icc_tbl  <- read_csv(file.path(paths$tables, "16_icc_summary.csv"),
                     show_col_types = FALSE) |>
  filter(!is.na(feature))

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
      p99_abs <= tau_99                                         ~ "unremarkable_p99",
      icc_player >= 0.35 & p99_abs > tau_99                    ~ "contamination-driven",
      icc_season >= 0.20 & !is.na(n_distinct_seasons) &
        n_distinct_seasons <= 2                                  ~ "season-driven",
      !is.na(recommended_transform) &
        recommended_transform != "none" & p99_abs > tau_99      ~ "skew-residual",
      TRUE                                                       ~ "mixed"
    )
  ) |>
  select(feature, group, p99_abs, tau_99_val, skewness, icc_player, icc_season,
         recommended_transform, n_distinct_seasons, verdict)

write_csv(tri, file.path(paths$tables, "29_residual_skew_triangulation.csv"))
message("Saved: 29_residual_skew_triangulation.csv")
message("\nVerdict counts:")
print(count(tri, verdict))

# ── Step 7 — Football context for top 20 worst cells ─────────────────────────

# season_index 1 = "2012/13", ..., 12 = "2023/24"
season_labels <- setNames(
  c("2012/13","2013/14","2014/15","2015/16","2016/17","2017/18",
    "2018/19","2019/20","2020/21","2021/22","2022/23","2023/24"),
  seq_len(S))

worst_top20 <- worst_tbl |>
  slice_head(n = 20) |>
  mutate(
    season_label = season_labels[as.character(season_index)],
    context = case_when(
      grepl("back_passes|^per90_passes$", feature) & abs(z_std) > 5 ~
        "Possible GK row missed by filter — verify player role",
      feature == "per90_smart_passes" & season_id == 8504 ~
        paste0("Reyes 2013/14 smart-passes: old model 14.4σ, new z=", round(z_std, 2)),
      feature == "per90_smart_passes" ~
        "High smart-pass season — creative playmaker or unusual role",
      grepl("goals|xg_shot|touch_in_box", feature) & z_std > 0 ~
        "Exceptional scoring/attacking season",
      grepl("interceptions|recoveries|clearances", feature) & z_std < -3 ~
        "Unusually low defensive output — injury, tactical shift, or role change",
      grepl("interceptions|recoveries", feature) & z_std > 3 ~
        "Exceptional defensive output season",
      grepl("dribbles|progressive_run|accelerations", feature) ~
        "Breakout dribbling/carrying season",
      grepl("crosses|shot_assists|key_passes", feature) & z_std > 3 ~
        "High-creativity wide/creative season",
      grepl("fouls|yellow", feature) ~
        "Disciplinary extreme — aggressive/fouled-heavily season",
      TRUE ~ "Manual review recommended"
    )
  )

context_lines <- c(
  "Stage 28 (K=3, P=48, N=4586, nu=4.902) — Top 20 worst cells by |z_std|",
  sprintf("tau_99 = %.4f (nu=%.4f). z|z| > tau_99 = model's own 99th-pct threshold.", tau_99, nu_hat),
  "",
  format(worst_top20 |> select(rank, player_name, season_label, feature, z_std, context),
         digits = 3)
)
writeLines(context_lines, file.path(paths$notes, "29_worst_cells_context.txt"))
message("Saved: 29_worst_cells_context.txt")

# ── Step 8 — Player factor scores and rankings ────────────────────────────────

ps <- read_csv(file.path(paths$tables,
                         "28_real_lowrank_a_diag_b_t_k3_posterior_summary.csv"),
               show_col_types = FALSE)

la_rows <- ps |>
  filter(grepl("^Lambda_a\\[", variable)) |>
  mutate(
    row_i = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\1", variable)),
    col_k = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\2", variable))
  )
L_mean <- matrix(0, P, K)
for (r in seq_len(nrow(la_rows))) L_mean[la_rows$row_i[r], la_rows$col_k[r]] <- la_rows$mean[r]

LLt      <- tcrossprod(L_mean)
eig      <- eigen(LLt, symmetric = TRUE)
V        <- eig$vectors[, 1:K]
var_exp  <- eig$values[1:K] / sum(eig$values[1:K])
message(sprintf("PCA of ΛΛ': PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%",
                100*var_exp[1], 100*var_exp[2], 100*var_exp[3]))

pc_scores <- A_mean %*% V    # I × K

player_info <- row_map |>
  select(player_index, player_id, player_name) |>
  distinct() |>
  arrange(player_index)

n_seasons_per_player <- row_map |>
  count(player_index, name = "n_seasons")

scores_tbl <- player_info |>
  left_join(n_seasons_per_player, by = "player_index") |>
  mutate(
    pc1_score = pc_scores[player_index, 1],
    pc2_score = pc_scores[player_index, 2],
    pc3_score = pc_scores[player_index, 3]
  )

write_csv(scores_tbl, file.path(paths$tables, "29_player_factor_scores.csv"))
message("Saved: 29_player_factor_scores.csv")

# Top / bottom 20 per PC
top_tbl <- bind_rows(lapply(1:K, function(k) {
  col_nm <- paste0("pc", k, "_score")
  bind_rows(
    scores_tbl |> slice_max(order_by = .data[[col_nm]], n = 20) |>
      mutate(pc = paste0("PC", k), direction = "top", rank_within = seq_len(n())),
    scores_tbl |> slice_min(order_by = .data[[col_nm]], n = 20) |>
      mutate(pc = paste0("PC", k), direction = "bottom",
             rank_within = seq_len(n()))
  )
}))
write_csv(top_tbl, file.path(paths$tables, "29_top_players_by_pc.csv"))
message("Saved: 29_top_players_by_pc.csv")

# Season factor trends from B_mean
season_pc_df <- data.frame(
  season_index = seq_len(S),
  season_pc1   = as.vector(B_mean %*% V[, 1]),
  season_pc2   = as.vector(B_mean %*% V[, 2]),
  season_pc3   = as.vector(B_mean %*% V[, 3])
) |>
  left_join(mo$season_map, by = "season_index") |>
  mutate(season_label = season_labels[as.character(season_index)])

write_csv(season_pc_df, file.path(paths$tables, "29_season_factor_trends.csv"))
message("Saved: 29_season_factor_trends.csv")

# ── Step 9 — Figures ──────────────────────────────────────────────────────────

save_plot <- function(p, fname, w = 10, h = 8) {
  ggplot2::ggsave(fname, p, width = w, height = h, dpi = 150)
  message("Saved: ", basename(fname))
}

verdict_colours <- c(
  "contamination-driven" = "#d73027",
  "season-driven"        = "#fc8d59",
  "skew-residual"        = "#fee090",
  "mixed"                = "#abd9e9",
  "unremarkable_p99"     = "#d0d0d0"
)

# Figure 1: residual p99 bar chart (horizontal)
p_bar <- tri |>
  mutate(feature = factor(feature, levels = rev(feat_smry$feature))) |>
  ggplot(aes(x = p99_abs, y = feature, fill = verdict)) +
  geom_col() +
  geom_vline(xintercept = tau_99, linetype = "dashed", colour = "black", linewidth = 0.8) +
  annotate("text", x = tau_99 + 0.02, y = 1, hjust = 0, size = 3,
           label = sprintf("tau99 = %.3f", tau_99)) +
  scale_fill_manual(values = verdict_colours, name = "Verdict") +
  labs(
    title   = "Stage 28 — Residual tail by feature (P=48, K=3, t-model)",
    subtitle = sprintf("Threshold tau99 = %.4f at nu = %.3f; dashed line = model's own 99th-pct", tau_99, nu_hat),
    x = "p99 of |z_std|", y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))

save_plot(p_bar, file.path(paths$figures, "29_residual_p99_barchart.png"), w = 12, h = 10)

# Figure 2: PC1 vs PC2 scatter
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
label_players <- scores_tbl |>
  bind_rows(
    scores_tbl |> slice_max(abs(pc1_score), n = 12),
    scores_tbl |> slice_max(abs(pc2_score), n = 12)
  ) |>
  distinct()

p_scatter <- ggplot(scores_tbl, aes(x = pc1_score, y = pc2_score)) +
  geom_point(aes(size = n_seasons), alpha = 0.35, colour = "#2166ac") +
  scale_size_continuous(range = c(0.5, 3), name = "Seasons") +
  labs(
    title    = sprintf("Stage 28 — Player archetypes: PC1 (%.1f%%) vs PC2 (%.1f%%)",
                       100*var_exp[1], 100*var_exp[2]),
    subtitle = "PC1: Offensive output (right) vs Defensive recovery (left)  |  PC2: Wide/creative (top)",
    x = sprintf("PC1 score (%.1f%% common variance)", 100*var_exp[1]),
    y = sprintf("PC2 score (%.1f%% common variance)", 100*var_exp[2])
  ) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  theme_minimal(base_size = 10)

if (has_ggrepel) {
  p_scatter <- p_scatter +
    ggrepel::geom_text_repel(data = label_players,
                              aes(label = player_name),
                              size = 2.5, max.overlaps = 25,
                              segment.colour = "grey60", segment.size = 0.3)
} else {
  thresh_label <- 2 * sd(scores_tbl$pc1_score)
  label_df <- scores_tbl |> filter(abs(pc1_score) > thresh_label | abs(pc2_score) > thresh_label)
  p_scatter <- p_scatter +
    geom_text(data = label_df, aes(label = player_name), size = 2.5, hjust = -0.1)
}

save_plot(p_scatter, file.path(paths$figures, "29_factor_score_scatter_pc1_pc2.png"),
          w = 11, h = 8)

# Figures 3a-c: top/bottom player tables per PC
make_pc_table_plot <- function(pc_num) {
  col_nm <- paste0("pc", pc_num, "_score")
  top15 <- scores_tbl |>
    slice_max(order_by = .data[[col_nm]], n = 15) |>
    arrange(desc(.data[[col_nm]])) |>
    mutate(rank_disp = seq_len(n()), x_pos = 0.05,
           y_pos = 15L - seq_len(n()) + 1L)
  bot15 <- scores_tbl |>
    slice_min(order_by = .data[[col_nm]], n = 15) |>
    arrange(.data[[col_nm]]) |>
    mutate(rank_disp = seq_len(n()), x_pos = 0.55,
           y_pos = 15L - seq_len(n()) + 1L)
  tbl_df <- bind_rows(top15, bot15) |>
    mutate(score_label = sprintf("%.3f", .data[[col_nm]]),
           label = paste0(rank_disp, ". ", player_name,
                          "  ", score_label, "  (", n_seasons, "s)"))

  ggplot(tbl_df, aes(x = x_pos, y = y_pos, label = label)) +
    geom_text(hjust = 0, size = 2.8) +
    annotate("text", x = 0.05, y = 16, hjust = 0, fontface = "bold", size = 3.5,
             label = paste0("Highest PC", pc_num)) +
    annotate("text", x = 0.55, y = 16, hjust = 0, fontface = "bold", size = 3.5,
             label = paste0("Lowest PC", pc_num)) +
    xlim(0, 1.1) + ylim(0, 17) +
    labs(title = sprintf("Stage 28 — PC%d archetype players (%.1f%% of shared variance)",
                         pc_num, 100 * var_exp[pc_num]),
         x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(plot.title  = element_text(hjust = 0.5, size = 11, face = "bold"),
          plot.margin = margin(8, 8, 8, 8))
}

for (k in 1:K) {
  p_tbl <- make_pc_table_plot(k)
  save_plot(p_tbl, file.path(paths$figures, sprintf("29_top_players_table_pc%d.png", k)),
            w = 10, h = 9)
}

# Figure 4: season factor trends
p_season <- season_pc_df |>
  pivot_longer(cols = starts_with("season_pc"), names_to = "pc", values_to = "value") |>
  mutate(pc = recode(pc,
    "season_pc1" = sprintf("PC1 (%.1f%%)", 100*var_exp[1]),
    "season_pc2" = sprintf("PC2 (%.1f%%)", 100*var_exp[2]),
    "season_pc3" = sprintf("PC3 (%.1f%%)", 100*var_exp[3])
  )) |>
  ggplot(aes(x = season_label, y = value, colour = pc, group = pc)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_brewer(palette = "Set1", name = "Factor") +
  labs(
    title    = "Stage 28 — Season-level factor trends (B_mean projected onto PCA axes)",
    subtitle = "Positive values: season shifted league toward that archetype direction",
    x = "Season", y = "Season effect on PC axis"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_season, file.path(paths$figures, "29_season_factor_trends.png"), w = 11, h = 5)

# ── Step 10 — Comparison with old Task A ─────────────────────────────────────

old_tau <- qt(0.995, df = 2.95) / sqrt(2.95 / (2.95 - 2))
p_ddw   <- match("rate_defensive_duels_won", feature_names)
frac_ddw_new <- mean(abs(z_std[, p_ddw]) > tau_99)

p_smp      <- match("per90_smart_passes", feature_names)
reyes_rows <- which(row_map$player_name %in% c("R. Reyes", "Reyes", "R.Reyes") &
                    row_map$season_id == 8504)
reyes_z    <- if (length(reyes_rows) > 0)
  max(abs(z_std[reyes_rows, p_smp])) else NA_real_

total_exceed_new <- sum(abs(z_std) > tau_99)
total_cells      <- as.double(N) * P

comparison_lines <- c(
  "=== Old Task A vs Stage 28 comparison ===",
  sprintf("  Old tau99 (nu=2.95): %.4f", old_tau),
  sprintf("  New tau99 (nu=%.4f): %.4f", nu_hat, tau_99),
  sprintf("  Lighter tails (higher nu) → lower threshold → model more surprised by extremes"),
  "",
  sprintf("  rate_defensive_duels_won frac > tau99: %.4f (old: GK-contaminated, should be lower)",
          frac_ddw_new),
  sprintf("  Reyes 2013/14 per90_smart_passes: %s (old: 14.4σ)",
          if (!is.na(reyes_z)) sprintf("%.2fσ", reyes_z) else "not found in rebuilt data"),
  sprintf("  Total (n,p) cells > tau99: %d / %.0f (%.2f%%)",
          total_exceed_new, total_cells, 100 * total_exceed_new / total_cells),
  sprintf("  Expected under t(nu): ~%.2f%%", 100 * 2 * pt(-tau_99, df = nu_hat)),
  ""
)

message(paste(comparison_lines, collapse = "\n"))
writeLines(comparison_lines, file.path(paths$notes, "29_vs_old_task_a.txt"))
message("Saved: 29_vs_old_task_a.txt")

# ── Final model comparison table ──────────────────────────────────────────────

loo_tbl <- read_csv(file.path(paths$tables, "C_loo_k2_vs_k3.csv"), show_col_types = FALSE)
elpd_k2 <- loo_tbl |> filter(model == "K2_stage18") |> pull(elpd_loo)
elpd_k3 <- loo_tbl |> filter(model == "K3_stage28") |> pull(elpd_loo)
se_k2   <- loo_tbl |> filter(model == "K2_stage18") |> pull(se_loo)
se_k3   <- loo_tbl |> filter(model == "K3_stage28") |> pull(se_loo)
d_elpd  <- loo_tbl |> filter(model == "K3_stage28") |> pull(delta_elpd)
se_d    <- loo_tbl |> filter(model == "K3_stage28") |> pull(se_delta)

comparison_tbl <- tibble(
  comparison            = c(
    "Low-rank K=2 vs diagonal K=0",
    "t-model vs Normal, K=2 (old data)",
    "Low-rank K=3 vs K=2 (Pathfinder)",
    "Low-rank K=3 vs K=2 (PSIS-LOO)",
    "Low-rank K=3 vs diagonal K=0"
  ),
  delta_elpd = c(72881, 3724, 12125, round(d_elpd, 1), NA_real_),
  se_delta   = c(NA, NA, NA, round(se_d, 1), NA_real_),
  n_se       = c(NA, ">20", NA, round(d_elpd / se_d, 2), NA_real_),
  source     = c(
    "K-fold Pathfinder (old P=53)",
    "PSIS-LOO old data",
    "Pathfinder K=0-4 (rebuilt P=48)",
    "PSIS-LOO rebuilt P=48, t-model",
    "Stage 20"
  ),
  notes = c(
    "Stage 24; Pareto-k unreliable — K-fold is primary",
    "Stage 18 vs stage 10; pre-rebuild N=4944/P=53",
    "Stage 17 rebuilt; Reversal at K=4 → K*=3",
    "C_loo_k2_vs_k3.csv; 99.4% Pareto-k>0.5 — directional only",
    "PENDING — Block E (Stage 03 diagonal refit still running)"
  )
)
write_csv(comparison_tbl, file.path(paths$tables, "final_model_comparison.csv"))
message("Saved: final_model_comparison.csv")

message("\nBlock A (stage 29) complete.")
