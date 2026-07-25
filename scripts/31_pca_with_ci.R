# Stage 31 — PCA with credible intervals on stage-28 posterior (Tasks 3+4, Session 11)
#
# Víctor: "Λ ni lo reportaría" — raw loadings are identified only up to rotation.
# This script replaces raw-Λ reporting with:
#   (a) PCA of posterior-mean ΛΛ' with per-entry q5/q95 CIs on factor loadings
#   (b) Same PCA of Σ_a = ΛΛ' + Ψ_a (adds per-feature uniqueness)
#   (c) CI-filtered top loaders: reliable = q5 and q95 have same sign
#
# Does NOT call as_cmdstan_fit() — reads from posterior summary CSV only.
# Outputs feed professor_summary.qmd (replace raw-Λ heatmap + top-loaders table).

source("src/bootstrap.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})
source("src/loading_visualization.R")

# ── 1. Load posterior summary ─────────────────────────────────────────────────

ps_path <- file.path(paths$tables, "28_real_lowrank_a_diag_b_t_k3_posterior_summary.csv")
if (!file.exists(ps_path)) stop("Posterior summary not found: ", ps_path)
ps <- read_csv(ps_path, show_col_types = FALSE)

mo           <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names          # length 48
K <- 3L; P <- length(feature_names)

# ── 2. Extract Lambda_a and psi_a ─────────────────────────────────────────────

la <- ps |>
  filter(grepl("^Lambda_a\\[", variable)) |>
  mutate(
    row_i = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\1", variable)),
    col_k = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\2", variable)),
    feature = feature_names[row_i]
  )

# Build L_mean (P×K) from posterior means
L_mean <- matrix(0, P, K)
for (r in seq_len(nrow(la))) L_mean[la$row_i[r], la$col_k[r]] <- la$mean[r]

psi_a <- ps |>
  filter(grepl("^psi_a\\[", variable)) |>
  arrange(as.integer(sub("psi_a\\[(\\d+)\\]", "\\1", variable))) |>
  pull(mean)
stopifnot(length(psi_a) == P)

cat(sprintf("Lambda_a entries: %d  (P=%d, K=%d)\n", nrow(la), P, K))
cat(sprintf("psi_a entries: %d\n", length(psi_a)))

# ── 3. PCA of ΛΛ' (common variance) ──────────────────────────────────────────
# SVD of L_mean: L_mean = U D V' where U (P×K) = left singular vectors,
# V (K×K) = right singular vectors, D (K×K) = diagonal of singular values.
# Then ΛΛ' = U D² U', so eigenvectors of ΛΛ' = U and eigenvalues = D².
# PCA loadings = U * D = L_mean * V  (P×K rotation from Lambda_a to PCA space).
# CI propagation: pca_lo = la_q5 %*% V  (P×K %*% K×K = P×K) — uses mean rotation.

LLt      <- tcrossprod(L_mean)             # P×P, used for Σ_a in section 5
svd_L    <- svd(L_mean, nu = K, nv = K)
U_pca    <- svd_L$u                        # P×K left singular vectors (= eigenvectors of ΛΛ')
V_right  <- svd_L$v                        # K×K right singular vectors
svals    <- svd_L$d                        # length K singular values
var_LLt  <- svals^2                        # eigenvalues of ΛΛ'
var_share_LLt <- var_LLt / sum(svals^2)               # fraction of total ΛΛ' trace variance

cat("\nΛΛ' PCA variance shares:\n")
for (k in 1:K)
  cat(sprintf("  PC%d: %.1f%%\n", k, 100 * var_share_LLt[k]))

# PCA loadings: U * D  (equivalent to L_mean %*% V_right)
loadings_LLt <- U_pca %*% diag(svals)    # P×K

# Fix sign: anchor each PC on its largest-magnitude feature (positive)
for (k in 1:K) {
  anchor <- which.max(abs(loadings_LLt[, k]))
  if (loadings_LLt[anchor, k] < 0) {
    loadings_LLt[, k] <- -loadings_LLt[, k]
    V_right[, k]      <- -V_right[, k]   # keep V_right consistent for CI propagation
  }
}

# ── 4. CI on PCA loadings via posterior q5/q95 of Lambda_a ───────────────────
# L_pca = L_mean %*% V_right  (exact linear map when V_right is held fixed at mean)
# So: pca_lo ≈ la_q5 %*% V_right  (P×K %*% K×K = P×K) — marginal, mean-rotation approx.

la_q5  <- matrix(0, P, K)
la_q95 <- matrix(0, P, K)
for (r in seq_len(nrow(la))) {
  la_q5[la$row_i[r],  la$col_k[r]] <- la$q5[r]
  la_q95[la$row_i[r], la$col_k[r]] <- la$q95[r]
}

# pmin/pmax across the two marginal extremes gives correct lo/hi after sign-mixing rotation
pca_lo <- pmin(la_q5 %*% V_right, la_q95 %*% V_right)  # P×K lower CI bound
pca_hi <- pmax(la_q5 %*% V_right, la_q95 %*% V_right)  # P×K upper CI bound

# ── 5. PCA of Σ_a = ΛΛ' + Ψ_a ────────────────────────────────────────────────

Sigma_a    <- LLt + diag(psi_a^2)
eg_Sa      <- eigen(Sigma_a, symmetric = TRUE)
V_Sa       <- eg_Sa$vectors[, 1:K]
var_Sa     <- eg_Sa$values[1:K]
var_share_Sa <- var_Sa / sum(eg_Sa$values)

cat("\nΣ_a PCA variance shares:\n")
for (k in 1:K)
  cat(sprintf("  PC%d: %.1f%%\n", k, 100 * var_share_Sa[k]))

loadings_Sa <- V_Sa %*% diag(sqrt(pmax(var_Sa, 0)))
for (k in 1:K) {
  anchor <- which.max(abs(loadings_Sa[, k]))
  if (loadings_Sa[anchor, k] < 0) loadings_Sa[, k] <- -loadings_Sa[, k]
}

# ── 6. CI-filtered top loaders table ─────────────────────────────────────────

groups <- as.data.frame(feature_group_lookup())

pca_ci_tbl <- do.call(rbind, lapply(1:K, function(k) {
  tibble(
    pc      = paste0("PC", k),
    feature = feature_names,
    group   = groups$group[match(feature_names, groups$feature)],
    loading = loadings_LLt[, k],
    ci_lo   = pca_lo[, k],
    ci_hi   = pca_hi[, k]
  ) |>
    mutate(reliable = (ci_lo > 0 & ci_hi > 0) | (ci_lo < 0 & ci_hi < 0))
}))

write_csv(pca_ci_tbl, file.path(paths$tables, "31_pca_loading_ci.csv"))
message("Saved: 31_pca_loading_ci.csv")

# ── 7. Σ_a comparison table ───────────────────────────────────────────────────

sigma_a_tbl <- do.call(rbind, lapply(1:K, function(k) {
  tibble(
    pc           = paste0("PC", k),
    feature      = feature_names,
    group        = groups$group[match(feature_names, groups$feature)],
    loading_LLt  = loadings_LLt[, k],
    loading_Sa   = loadings_Sa[, k]
  )
}))

eigenval_tbl <- tibble(
  pc            = paste0("PC", 1:K),
  eigenval_LLt  = var_LLt,
  var_share_LLt = var_share_LLt,
  eigenval_Sa   = var_Sa,
  var_share_Sa  = var_share_Sa
)

write_csv(sigma_a_tbl,  file.path(paths$tables, "31_sigma_a_pca_loadings.csv"))
write_csv(eigenval_tbl, file.path(paths$tables, "31_sigma_a_pca.csv"))
message("Saved: 31_sigma_a_pca.csv, 31_sigma_a_pca_loadings.csv")

# ── 8. CI errorbar loading plot (reliable loaders highlighted) ────────────────

plot_ci <- pca_ci_tbl |>
  group_by(pc) |>
  slice_max(order_by = abs(loading), n = 15) |>
  ungroup() |>
  mutate(
    feature_short = gsub("_", " ", gsub("^per90_|^rate_", "", feature)),
    feature_f     = reorder(interaction(pc, feature_short), abs(loading)),
    alpha_val     = ifelse(reliable, 1.0, 0.35)
  )

p_ci <- ggplot(plot_ci, aes(x = loading, y = feature_f, colour = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, alpha = alpha_val),
                 height = 0.3, linewidth = 0.7) +
  geom_point(aes(alpha = alpha_val), size = 2) +
  scale_alpha_identity() +
  scale_colour_manual(values = group_colours(), name = "Feature group", na.value = "grey50") +
  facet_wrap(~pc, scales = "free_y", ncol = 3) +
  labs(
    title    = "PCA factor loadings with 90% credible intervals — stage 28 (K=3)",
    subtitle = paste0(
      "Solid = CI excludes zero (reliable loader).  ",
      sprintf("ΛΛ' variance shares: PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%",
              100*var_share_LLt[1], 100*var_share_LLt[2], 100*var_share_LLt[3])
    ),
    x = "PCA loading  (eigenvector × √eigenvalue)", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold"),
    axis.text.y     = element_text(size = 7)
  )

ggsave(file.path(paths$figures, "31_pca_loading_ci_plot.png"),
       p_ci, width = 13, height = 9, dpi = 150)
message("Saved: 31_pca_loading_ci_plot.png")

# ── 9. Eigenvalue share comparison: ΛΛ' vs Σ_a ───────────────────────────────

eig_long <- bind_rows(
  eigenval_tbl |> transmute(pc, var_share = var_share_LLt, Source = "ΛΛ' (common)"),
  eigenval_tbl |> transmute(pc, var_share = var_share_Sa,  Source = "Σ_a = ΛΛ' + Ψ_a (total player)")
)

p_eig <- ggplot(eig_long, aes(x = pc, y = 100 * var_share, fill = Source)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * var_share)),
            position = position_dodge(width = 0.6), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("ΛΛ' (common)" = "#2166ac",
                                "Σ_a = ΛΛ' + Ψ_a (total player)" = "#d73027"),
                    name = NULL) +
  expand_limits(y = max(100 * eig_long$var_share) + 5) +
  labs(
    title    = "Variance shares: ΛΛ' vs Σ_a per PC (stage 28, K=3)",
    subtitle = "Σ_a adds per-feature uniqueness Ψ_a. Shares are fraction of total variance of each matrix.",
    x = NULL, y = "Variance share (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

ggsave(file.path(paths$figures, "31_sigma_a_vs_llambda_comparison.png"),
       p_eig, width = 7, height = 5, dpi = 150)
message("Saved: 31_sigma_a_vs_llambda_comparison.png")

# ── 10. Print top reliable loaders per PC ────────────────────────────────────

cat("\n=== CI-filtered top loaders (reliable = CI excludes zero) ===\n")
for (k in 1:K) {
  cat(sprintf("\nPC%d (ΛΛ' share: %.1f%%):\n", k, 100 * var_share_LLt[k]))
  top <- pca_ci_tbl |>
    filter(pc == paste0("PC", k), reliable) |>
    arrange(desc(abs(loading))) |>
    head(8) |>
    mutate(feature_short = gsub("_", " ", gsub("^per90_|^rate_", "", feature)))
  print(top[, c("feature_short", "group", "loading", "ci_lo", "ci_hi")])
}

message("\nStage 31 complete.")
