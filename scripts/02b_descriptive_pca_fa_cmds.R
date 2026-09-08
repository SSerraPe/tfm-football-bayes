# Descriptive low-dimensional structure analysis (PCA / FA / CMDS) on the production
# modelling dataset (N=4586, P=48, Z-scaled + transformed in stage 02).
#
# This reimplements, on the corrected/production dataset, the exploratory analysis
# originally prototyped in tfm_playground/similarity_measures.Rmd (block-weighted
# distance construction) and tfm_playground/playground_BMDS.Rmd (PCA, parallel-analysis
# FA with oblimin rotation, classical MDS). See archive/exploratory_p23_similarity_analysis.R
# for the superseded p=23 version this replaces as the source for thesis Chapter 3
# (docs are tfm_latex/data_source.tex and tfm_latex/background.tex).
#
# Unlike the playground notebooks, this script does NOT re-apply log1p to the input
# matrix: stage 02 (scripts/02_prepare_model_objects.R) already applies feature-specific
# log1p/sqrt transforms *before* Z-scaling, so `model_objects$Y` is already on an
# appropriate scale (and, being Z-scored, contains negative values that a second log1p
# pass would turn into NaN). Block-weighting is applied on top of the already-scaled Y,
# exactly as in the original notebook's `col_w` step.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/02b_descriptive_pca_fa_cmds.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c("dplyr", "tidyr", "tibble", "ggplot2", "psych", "ggrepel", "scales"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

set.seed(pipeline$seed)

model_objects_path <- file.path(paths$processed, "model_objects.rds")
if (!file.exists(model_objects_path)) {
  stop("Missing model_objects.rds. Run stage 02 first.", call. = FALSE)
}
mo <- readRDS(model_objects_path)

Y <- mo$Y
meta <- mo$metadata |>
  mutate(season_id = as.factor(season_id))
stopifnot(nrow(Y) == nrow(meta))

N <- nrow(Y)
P <- ncol(Y)
feat_names <- colnames(Y)

# ---------------------------------------------------------------------------
# 1) Football blocks (same 5-block scheme as the original exploratory notebook)
# ---------------------------------------------------------------------------
blocks <- list(
  offense = c(
    "per90_goals", "per90_assists", "per90_xg_shot", "per90_xg_assist",
    "per90_shots", "per90_head_shots", "per90_touch_in_box"
  ),
  progression_creation = c(
    "per90_received_pass", "per90_passes", "per90_forward_passes", "per90_back_passes",
    "per90_lateral_passes", "per90_long_passes", "per90_smart_passes",
    "per90_progressive_passes", "per90_passes_to_final_third", "per90_key_passes",
    "per90_through_passes", "per90_crosses", "per90_shot_assists",
    "rate_successful_passes", "rate_successful_progressive_passes"
  ),
  carrying_1v1 = c(
    "per90_dribbles", "per90_progressive_run", "per90_accelerations",
    "rate_successful_dribbles"
  ),
  defense_press = c(
    "per90_losses", "per90_own_half_losses", "per90_dangerous_own_half_losses",
    "per90_recoveries", "per90_opponent_half_recoveries",
    "per90_dangerous_opponent_half_recoveries", "per90_interceptions",
    "per90_clearances", "per90_shots_blocked", "per90_sliding_tackles",
    "per90_dribbles_against", "per90_pressing_duels"
  ),
  # Discipline features (fouls, cards) are folded into duels_phys: the original
  # 5-block scheme has no discipline bucket, and fouls/cards are physical-
  # confrontation-adjacent.
  duels_phys = c(
    "per90_defensive_duels", "per90_offensive_duels", "per90_aerial_duels",
    "per90_loose_ball_duels", "rate_defensive_duels_won", "rate_offensive_duels_won",
    "rate_aerial_duels_won", "per90_fouls", "per90_fouls_suffered", "per90_yellow_cards"
  )
)

all_block_feats <- unlist(blocks, use.names = FALSE)
stopifnot(setequal(all_block_feats, feat_names))
stopifnot(!anyDuplicated(all_block_feats))

feat_to_block <- tibble(
  feature = all_block_feats,
  block = rep(names(blocks), times = vapply(blocks, length, integer(1)))
)

block_weights <- setNames(rep(1.0, length(blocks)), names(blocks))
col_w <- setNames(rep(1, P), feat_names)
for (b in names(blocks)) {
  feats_b <- blocks[[b]]
  col_w[feats_b] <- block_weights[[b]] / sqrt(length(feats_b))
}

X_w <- sweep(Y, 2, col_w[feat_names], `*`)

# ---------------------------------------------------------------------------
# 2) Distance matrices (Euclidean + cosine) on the block-weighted matrix
# ---------------------------------------------------------------------------
cosine_dist <- function(M) {
  M <- as.matrix(M)
  M[is.na(M)] <- 0
  norms <- sqrt(rowSums(M^2))
  norms[norms == 0] <- 1
  Mnorm <- M / norms
  S <- Mnorm %*% t(Mnorm)
  D <- 1 - S
  diag(D) <- 0
  as.dist(D)
}

D_euclid <- dist(X_w, method = "euclidean")
D_cosine <- cosine_dist(X_w)

# ---------------------------------------------------------------------------
# 3) Euclidean-vs-cosine agreement figure (N=900 stratified-by-season subset,
#    same construction as tfm_playground/playground_BMDS.Rmd's BMDS subset, but
#    only used here for a classical-methods comparison figure)
# ---------------------------------------------------------------------------
min_minutes_subset <- 900
N_target <- 900
meta0 <- meta |> mutate(row_id = row_number())
meta_f0 <- meta0 |> filter(minutes >= min_minutes_subset)
n_seasons <- n_distinct(meta_f0$season_id)
n_per_season <- ceiling(N_target / n_seasons)

idx_keep <- meta_f0 |>
  group_by(season_id) |>
  group_modify(~ slice_sample(.x, n = min(nrow(.x), n_per_season))) |>
  ungroup()
keep_rows <- sort(idx_keep$row_id)

D_euc_sub <- as.matrix(D_euclid)[keep_rows, keep_rows]
D_cos_sub <- as.matrix(D_cosine)[keep_rows, keep_rows]
d_euc_vec <- D_euc_sub[upper.tri(D_euc_sub)]
d_cos_vec <- D_cos_sub[upper.tri(D_cos_sub)]

cor_pearson_ec <- cor(d_euc_vec, d_cos_vec, method = "pearson")
cor_spearman_ec <- cor(d_euc_vec, d_cos_vec, method = "spearman")

n_pairs_avail <- length(d_euc_vec)
idx_sample <- sample.int(n_pairs_avail, min(30000, n_pairs_avail))
dist_compare_df <- tibble(
  d_euclid = d_euc_vec[idx_sample],
  d_coseno = d_cos_vec[idx_sample]
)

p_scatter_distances <- ggplot(dist_compare_df, aes(x = d_euclid, y = d_coseno)) +
  geom_point(alpha = 0.08, size = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 0.8) +
  theme_minimal() +
  labs(x = "Euclidean distance", y = "Cosine dissimilarity")

ggsave(file.path(paths$scripts, "..", "tfm_latex", "images", "scatter_distances.png"),
       p_scatter_distances, width = 6.5, height = 5, dpi = 300)

writeLines(
  c(
    sprintf("N (subset) = %d", length(keep_rows)),
    sprintf("n_pairs_sampled = %d", length(idx_sample)),
    sprintf("Pearson(d_euclid, d_cosine) = %.4f", cor_pearson_ec),
    sprintf("Spearman(d_euclid, d_cosine) = %.4f", cor_spearman_ec)
  ),
  file.path(paths$notes, "02b_distance_agreement.txt")
)

# ---------------------------------------------------------------------------
# 4) PCA on the block-weighted, already-scaled matrix
# ---------------------------------------------------------------------------
pca <- prcomp(X_w, center = FALSE, scale. = FALSE)
eig_vals <- pca$sdev^2
prop_var <- eig_vals / sum(eig_vals)
cum_var <- cumsum(prop_var)
q90 <- which(cum_var >= 0.90)[1]

pca_var_df <- tibble(PC = seq_along(prop_var), prop_var = prop_var, cum_var = cum_var)
write_csv_safe <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, path)
}
write_csv_safe(pca_var_df, file.path(paths$tables, "02b_pca_variance_explained.csv"))

p_scree_pca <- ggplot(pca_var_df |> dplyr::slice(1:min(30, P)), aes(x = PC, y = prop_var)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  geom_vline(xintercept = q90, linetype = "dashed", color = "firebrick", linewidth = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal() +
  labs(x = "Principal component", y = "Proportion of variance explained")
ggsave(file.path(paths$scripts, "..", "tfm_latex", "images", "scree_pca.png"),
       p_scree_pca, width = 7, height = 5, dpi = 300)

# PC1 vs PC2-5, coloured by minutes
scores <- as.data.frame(pca$x[, 1:5, drop = FALSE])
colnames(scores) <- paste0("PC", 1:5)
pca_df <- meta |> bind_cols(scores)

p_scatter_pca <- pca_df |>
  pivot_longer(cols = c(PC2, PC3, PC4, PC5), names_to = "component", values_to = "score") |>
  ggplot(aes(x = PC1, y = score, color = minutes)) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  facet_wrap(~component, scales = "free_y", nrow = 2) +
  theme_minimal() +
  labs(x = paste0("PC1 (", round(100 * prop_var[1], 1), "%)"), y = "Score", color = "Minutes played")
ggsave(file.path(paths$scripts, "..", "tfm_latex", "images", "scatter_pca.png"),
       p_scatter_pca, width = 7, height = 5, dpi = 300)

pc_min_cor <- pca_df |>
  summarise(
    cor_PC1_minutes = cor(PC1, minutes, use = "complete.obs"),
    cor_PC2_minutes = cor(PC2, minutes, use = "complete.obs"),
    cor_PC3_minutes = cor(PC3, minutes, use = "complete.obs")
  )
write_csv_safe(pc_min_cor, file.path(paths$tables, "02b_pca_minutes_correlation.csv"))

# Top-5 absolute loadings for the first 5 PCs
loadings_mat <- pca$rotation
N_PCS_TO_SHOW <- 5
TOP_N_LOADINGS <- 5
pca_loadings_long <- as.data.frame(loadings_mat[, 1:N_PCS_TO_SHOW, drop = FALSE]) |>
  rownames_to_column("feature") |>
  pivot_longer(cols = -feature, names_to = "PC", values_to = "loading") |>
  mutate(abs_loading = abs(loading))

top_pca_loadings <- pca_loadings_long |>
  group_by(PC) |>
  slice_max(order_by = abs_loading, n = TOP_N_LOADINGS, with_ties = FALSE) |>
  arrange(PC, desc(abs_loading)) |>
  ungroup()
write_csv_safe(top_pca_loadings, file.path(paths$tables, "02b_pca_top_loadings.csv"))

# ---------------------------------------------------------------------------
# 5) Factor analysis: parallel analysis + oblimin rotation
# ---------------------------------------------------------------------------
sds <- apply(X_w, 2, sd)
zero_sd_cols <- which(is.na(sds) | sds == 0)
X_fa <- if (length(zero_sd_cols) > 0) X_w[, -zero_sd_cols, drop = FALSE] else X_w

pa <- psych::fa.parallel(X_fa, fa = "fa", n.iter = 200, plot = FALSE)
n_factors <- pa$nfact
pa_df <- tibble(factor = seq_along(pa$fa.values), observed = pa$fa.values, simulated = pa$fa.sim)
write_csv_safe(pa_df, file.path(paths$tables, "02b_fa_parallel_analysis.csv"))

p_scree_fa <- ggplot(pa_df |> dplyr::slice(1:min(20, nrow(pa_df))), aes(x = factor)) +
  geom_line(aes(y = observed, color = "Observed data"), linewidth = 0.8) +
  geom_point(aes(y = observed, color = "Observed data"), size = 2) +
  geom_line(aes(y = simulated, color = "Simulated (95th pct)"), linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(y = simulated, color = "Simulated (95th pct)"), size = 2) +
  geom_vline(xintercept = n_factors, linetype = "dashed", color = "firebrick", linewidth = 0.6) +
  scale_color_manual(values = c("Observed data" = "steelblue", "Simulated (95th pct)" = "gray40")) +
  theme_minimal() +
  labs(x = "Factor", y = "Eigenvalue", color = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(paths$scripts, "..", "tfm_latex", "images", "scree_fa.png"),
       p_scree_fa, width = 7, height = 5, dpi = 300)

fit_fa <- function(k) psych::fa(X_fa, nfactors = k, rotate = "oblimin", fm = "ml", scores = "regression")
k_compare_lo <- max(2, n_factors - 2)
fa_lo <- fit_fa(k_compare_lo)
fa_hi <- fit_fa(n_factors)

fit_tbl <- tibble(
  k = c(k_compare_lo, n_factors),
  TLI = c(fa_lo$TLI, fa_hi$TLI),
  RMSEA = c(fa_lo$RMSEA[1], fa_hi$RMSEA[1]),
  RMSR = c(fa_lo$rms, fa_hi$rms),
  BIC = c(fa_lo$BIC, fa_hi$BIC)
)
write_csv_safe(fit_tbl, file.path(paths$tables, "02b_fa_fit_comparison.csv"))

L <- unclass(fa_hi$loadings)
storage.mode(L) <- "double"
load_long <- as.data.frame(L) |>
  rownames_to_column("feature") |>
  pivot_longer(-feature, names_to = "factor", values_to = "loading") |>
  mutate(abs_loading = abs(loading)) |>
  filter(abs_loading >= 0.30) |>
  arrange(factor, desc(abs_loading))
write_csv_safe(load_long, file.path(paths$tables, "02b_fa_loadings.csv"))

fa_cor_tbl <- as.data.frame(fa_hi$Phi)
write_csv_safe(rownames_to_column(fa_cor_tbl, "factor"), file.path(paths$tables, "02b_fa_factor_correlations.csv"))

# Panels are ordered by variance explained (fa_hi$Vaccounted["Proportion Var", ]) and
# labeled with the short interpretation given in the surrounding prose (background.tex,
# "Ordered by the variance each factor accounts for..."). If the FA is refit and the
# ML-column-to-archetype correspondence shifts, this vector must be re-checked against
# fa_hi$Vaccounted and updated to match.
factor_var_order <- names(sort(fa_hi$Vaccounted["Proportion Var", ], decreasing = TRUE))
factor_interpretation <- c(
  "chance creation & assists", "goal threat & finishing", "carrying & offensive duelling",
  "passing volume & circulation", "defensive duel-winning", "defensive duel exposure",
  "progressive & vertical passing", "physical engagement & possession risk",
  "high/pressing recoveries", "passing accuracy vs. crossing risk", "deep build-up circulation"
)
wrap_label <- function(x, width = 20) paste(strwrap(x, width = width), collapse = "\n")
factor_display <- setNames(
  paste0("F", seq_along(factor_var_order), ": ", vapply(factor_interpretation, wrap_label, character(1))),
  factor_var_order
)

load_long <- load_long |>
  mutate(
    factor_label = factor(factor_display[factor], levels = unname(factor_display[factor_var_order])),
    feature_label = gsub("per90_|rate_", "", feature),
    feature_label = gsub("_", " ", feature_label)
  )

p_loadings_fa <- ggplot(load_long, aes(x = reorder(feature_label, abs_loading), y = loading,
                                        fill = ifelse(loading > 0, "Positive", "Negative"))) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray30") +
  coord_flip() +
  facet_wrap(~factor_label, scales = "free_y", ncol = 4) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3), labels = scales::label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "firebrick")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    panel.spacing = unit(1, "lines")
  ) +
  labs(x = NULL, y = "Loading", fill = NULL)
ggsave(file.path(paths$scripts, "..", "tfm_latex", "images", "loadings_fa.png"),
       p_loadings_fa, width = 15, height = 10, dpi = 300)

# ---------------------------------------------------------------------------
# 6) Classical MDS (Euclidean + cosine), full sample; extremes from the SAME
#    (full-sample) embedding used for the plotted points — the original notebook
#    computed extremes from a *different* subset-only embedding and overlaid them
#    on the full-sample scatter, a coordinate-system mismatch this version fixes.
# ---------------------------------------------------------------------------
mds_cosine <- cmdscale(as.matrix(D_cosine), k = 2, eig = TRUE)
coords_mds <- as.data.frame(mds_cosine$points)
colnames(coords_mds) <- c("Dim1", "Dim2")
mds_df <- meta |> bind_cols(coords_mds)

p_cmds_cosine <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = minutes)) +
  geom_point(alpha = 0.4, size = 1.0) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(x = "Dimension 1", y = "Dimension 2", color = "Minutes played")
ggsave(file.path(paths$scripts, "..", "tfm_latex", "images", "cmds_cosine.png"),
       p_cmds_cosine, width = 7, height = 5.5, dpi = 300)

k_ext <- 5
top_bottom <- function(df, var, k) {
  bind_rows(
    df |> arrange(desc(.data[[var]])) |> slice(1:k) |> mutate(extreme = paste0("High ", var)),
    df |> arrange(.data[[var]]) |> slice(1:k) |> mutate(extreme = paste0("Low ", var))
  )
}
extremes <- bind_rows(top_bottom(mds_df, "Dim1", k_ext), top_bottom(mds_df, "Dim2", k_ext)) |>
  mutate(tag = paste0(player_name, "\n(", season_id, ")"))
write_csv_safe(
  extremes |> select(player_name, season_id, minutes, Dim1, Dim2, extreme),
  file.path(paths$tables, "02b_cmds_extreme_observations.csv")
)

p_cmds_extremes <- mds_df |>
  ggplot(aes(x = Dim1, y = Dim2)) +
  geom_point(alpha = 0.25, size = 0.8, color = "gray60") +
  geom_point(data = extremes, aes(color = extreme), size = 2.0, alpha = 0.9) +
  ggrepel::geom_text_repel(data = extremes, aes(label = tag, color = extreme),
                            size = 2.5, max.overlaps = 20, show.legend = FALSE) +
  theme_minimal() +
  labs(x = "Dimension 1", y = "Dimension 2", color = "Extreme")
ggsave(file.path(paths$scripts, "..", "tfm_latex", "images", "cmds_extreme_obs.png"),
       p_cmds_extremes, width = 7.5, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# 7) Summary note
# ---------------------------------------------------------------------------
writeLines(
  c(
    sprintf("N = %d, P = %d", N, P),
    sprintf("PCA: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%", 100 * prop_var[1], 100 * prop_var[2], 100 * prop_var[3]),
    sprintf("PCA cumulative: PC1-2=%.1f%%, PC1-5=%.1f%%, PC1-10=%.1f%%",
            100 * cum_var[2], 100 * cum_var[5], 100 * cum_var[10]),
    sprintf("PCs needed for 90%% cumulative variance: %d", q90),
    sprintf("FA parallel analysis suggested factors: %d", n_factors),
    sprintf("FA fit comparison: k=%d (TLI=%.3f, RMSEA=%.3f, RMSR=%.4f, BIC=%.1f) vs k=%d (TLI=%.3f, RMSEA=%.3f, RMSR=%.4f, BIC=%.1f)",
            k_compare_lo, fa_lo$TLI, fa_lo$RMSEA[1], fa_lo$rms, fa_lo$BIC,
            n_factors, fa_hi$TLI, fa_hi$RMSEA[1], fa_hi$rms, fa_hi$BIC),
    sprintf("Distance agreement (N=%d subset): Pearson=%.4f, Spearman=%.4f",
            length(keep_rows), cor_pearson_ec, cor_spearman_ec)
  ),
  file.path(paths$notes, "02b_descriptive_summary.txt")
)

message("Stage 02b complete. See outputs/notes/02b_descriptive_summary.txt for headline numbers.")
