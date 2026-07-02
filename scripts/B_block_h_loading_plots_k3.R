# Block H support: generate K=3 loading plots from stage 28 posterior summary.
# Produces:
#   28_lambda_a_heatmap.png   — raw LT loadings (P=48 × K=3)
#   28_pca_loading_heatmap.png — PCA-rotated loadings
#   28_factor_group_radar.png  — radar by feature group

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/B_block_h_loading_plots_k3.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "ggplot2", "tidyr"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(tidyr); library(tibble); library(forcats)
})

model_objects <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- model_objects$variable_names
P             <- length(feature_names)
K             <- 3L

ps_path <- file.path(paths$tables, "28_real_lowrank_a_diag_b_t_k3_posterior_summary.csv")
ps <- read_csv(ps_path, show_col_types = FALSE)

groups <- feature_group_lookup()

save_plot <- function(p, fname, w = 12, h = 14) {
  ggplot2::ggsave(fname, p, width = w, height = h, dpi = 150)
  message("Saved: ", basename(fname))
}

# ── Extract Lambda_a posterior mean matrix ───────────────────────────────────

la_rows <- ps[grepl("^Lambda_a\\[", ps$variable), ]
la_rows <- la_rows %>%
  mutate(
    row_i = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\1", variable)),
    col_k = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\2", variable))
  )

L_mean <- matrix(0, P, K)
for (r in seq_len(nrow(la_rows))) {
  L_mean[la_rows$row_i[r], la_rows$col_k[r]] <- la_rows$mean[r]
}
colnames(L_mean) <- paste0("Factor ", 1:K)
rownames(L_mean) <- feature_names

message("Lambda_a mean matrix built: ", P, " × ", K)

# ── 1. Raw Lambda_a heatmap ───────────────────────────────────────────────────

tbl_raw <- tibble(
  feature = rep(feature_names, K),
  factor  = rep(paste0("Factor ", 1:K), each = P),
  loading = as.vector(L_mean)
) %>%
  left_join(groups, by = "feature") %>%
  mutate(
    group   = coalesce(group, "Other"),
    group   = factor(group, levels = group_order()),
    feature = factor(feature, levels = feature_names),
    loading_lt = if_else(
      as.integer(feature) >= as.integer(factor(factor, levels = paste0("Factor ", 1:K))),
      loading, NA_real_
    )
  )

# For LT structure: only show Lambda_a[i,k] where i >= k
tbl_lt <- tbl_raw %>%
  mutate(
    k_idx = as.integer(sub("Factor ", "", factor)),
    i_idx = as.integer(feature),
    visible = i_idx >= k_idx
  ) %>%
  mutate(loading_plot = if_else(visible, loading, NA_real_))

p_heatmap <- ggplot(tbl_lt, aes(x = factor, y = fct_rev(feature), fill = loading_plot)) +
  geom_tile(colour = "grey90", linewidth = 0.2) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4",
                       midpoint = 0, na.value = "grey95",
                       name = "Loading\n(posterior mean)") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Stage 28 — Factor loadings Λ (K=3, P=48, t-model)",
    subtitle = "Lower-triangular with positive diagonal; upper triangle fixed at 0",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7, face = "bold"),
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_text(size = 8),
    legend.position = "bottom"
  )

save_plot(p_heatmap, file.path(paths$figures, "28_lambda_a_heatmap.png"), w = 9, h = 16)

# ── 2. PCA of ΛΛ' for canonical rotation ─────────────────────────────────────

LLt <- tcrossprod(L_mean)
eig <- eigen(LLt, symmetric = TRUE)
n_pc <- K
V <- eig$vectors[, 1:n_pc, drop = FALSE]   # P × K PCA loadings

var_explained <- eig$values[1:n_pc] / sum(eig$values[1:n_pc])
message("PCA of ΛΛ' — variance explained by each PC:")
for (k in seq_len(n_pc)) message(sprintf("  PC%d: %.1f%%", k, 100 * var_explained[k]))

tbl_pca <- tibble(
  feature = rep(feature_names, n_pc),
  pc      = rep(paste0("PC", 1:n_pc, " (", round(100 * var_explained, 1), "%)"), each = P),
  loading = as.vector(V)
) %>%
  left_join(groups, by = "feature") %>%
  mutate(
    group   = coalesce(group, "Other"),
    group   = factor(group, levels = group_order()),
    feature = factor(feature, levels = feature_names)
  )

p_pca_heatmap <- ggplot(tbl_pca, aes(x = pc, y = fct_rev(feature), fill = loading)) +
  geom_tile(colour = "grey90", linewidth = 0.2) +
  scale_fill_gradient2(low = "#d73027", mid = "white", high = "#4575b4",
                       midpoint = 0, name = "PCA loading") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = "Stage 28 — PCA of ΛΛ' (canonical factor rotation, K=3)",
    subtitle = "Eigenvectors of posterior-mean ΛΛ' give rotation-invariant factor interpretation",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, size = 7, face = "bold"),
    axis.text.y  = element_text(size = 6),
    axis.text.x  = element_text(size = 7),
    legend.position = "bottom"
  )

save_plot(p_pca_heatmap, file.path(paths$figures, "28_pca_loading_heatmap.png"), w = 9, h = 16)

# ── 3. Group-level radar (mean absolute loading per group per PC) ─────────────

if (requireNamespace("fmsb", quietly = TRUE)) {
  radar_df <- tbl_pca %>%
    group_by(pc, group) %>%
    summarise(mean_abs = mean(abs(loading)), .groups = "drop")

  p_radar <- ggplot(radar_df, aes(x = group, y = mean_abs, colour = pc, group = pc)) +
    geom_line() + geom_point(size = 2.5) +
    scale_colour_brewer(palette = "Set1", name = "Factor") +
    coord_polar() +
    labs(
      title = "Stage 28 — Mean absolute PCA loading by feature group",
      x = NULL, y = "Mean |loading|"
    ) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(size = 8))
  save_plot(p_radar, file.path(paths$figures, "28_factor_group_radar.png"), w = 8, h = 8)
} else {
  # Bar chart alternative
  radar_df <- tbl_pca %>%
    group_by(pc, group) %>%
    summarise(mean_abs = mean(abs(loading)), .groups = "drop")
  p_bar <- ggplot(radar_df, aes(x = group, y = mean_abs, fill = pc)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Set1", name = "PC") +
    labs(title = "Stage 28 — Mean |PCA loading| by feature group", x = NULL, y = "Mean |loading|") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  save_plot(p_bar, file.path(paths$figures, "28_factor_group_radar.png"), w = 10, h = 6)
}

message("Block H loading plots complete.")
