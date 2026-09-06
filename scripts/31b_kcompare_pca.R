# Stage 31b — Per-K PCA with credible intervals (Session 14).
#
# K-parameterised variant of scripts/31_pca_with_ci.R. Reads a posterior
# summary CSV, computes:
#   * PCA of posterior-mean Lambda_a Lambda_a' (common variance shares)
#   * PCA of Sigma_a = Lambda_a Lambda_a' + Psi_a  (total-player shares)
#   * 90% CIs on rotated loadings via linear map at posterior-mean rotation
#     (documented mean-rotation approximation, exactly as stage 31)
#   * Reliable-loader count per PC (CI excludes zero)
#
# Selects the fit via BFA_K env var (2, 3, or 4). Maps:
#   K=2 -> 18_real_lowrank_a_diag_b_t
#   K=3 -> 28_real_lowrank_a_diag_b_t_k3
#   K=4 -> 21_real_k4_t
#
# Writes:
#   outputs/tables/31b_k{K}_variance_shares.csv
#   outputs/tables/31b_k{K}_reliable_counts.csv
#   outputs/tables/31b_k{K}_pca_loading_ci.csv
#   outputs/tables/31b_kcompare_variance_shares.csv (appended across K when all K done)
#   outputs/tables/31b_kcompare_reliable_counts.csv (appended across K when all K done)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/31b_kcompare_pca.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))
source("src/loading_visualization.R")

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
  stop("Posterior summary not found: ", ps_path,
       "\n  (For K=", K, ", fit_id=", fit_id, ". ",
       "Check that stage ", substr(fit_id, 1, 2), " has completed.)")

message(sprintf("Stage 31b: K=%d, fit_id=%s", K, fit_id))

mo <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names
P <- length(feature_names)

# ── Load posterior summary ────────────────────────────────────────────────────

ps <- read_csv(ps_path, show_col_types = FALSE)

la <- ps |>
  filter(grepl("^Lambda_a\\[", variable)) |>
  mutate(
    row_i = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\1", variable)),
    col_k = as.integer(sub("Lambda_a\\[(\\d+),(\\d+)\\]", "\\2", variable))
  )
stopifnot(nrow(la) == P * K)

# Build L_mean, L_q5, L_q95
L_mean <- matrix(0, P, K)
L_q5   <- matrix(0, P, K)
L_q95  <- matrix(0, P, K)
for (r in seq_len(nrow(la))) {
  L_mean[la$row_i[r], la$col_k[r]] <- la$mean[r]
  L_q5  [la$row_i[r], la$col_k[r]] <- la$q5[r]
  L_q95 [la$row_i[r], la$col_k[r]] <- la$q95[r]
}

psi_a <- ps |>
  filter(grepl("^psi_a\\[", variable)) |>
  arrange(as.integer(sub("psi_a\\[(\\d+)\\]", "\\1", variable))) |>
  pull(mean)
stopifnot(length(psi_a) == P)

cat(sprintf("Lambda_a: %d entries (P=%d, K=%d), psi_a: %d entries\n",
            nrow(la), P, K, length(psi_a)))

# ── PCA of Lambda_a Lambda_a' (common) ────────────────────────────────────────

svd_L         <- svd(L_mean, nu = K, nv = K)
U_pca         <- svd_L$u                        # P x K
V_right       <- svd_L$v                        # K x K
svals         <- svd_L$d                        # length K
var_LLt       <- svals^2                        # eigenvalues
var_share_LLt <- var_LLt / sum(var_LLt)

# Sign-fix: each PC's max-magnitude feature is positive
loadings_LLt <- U_pca %*% diag(svals, K)
for (k in 1:K) {
  anchor <- which.max(abs(loadings_LLt[, k]))
  if (loadings_LLt[anchor, k] < 0) {
    loadings_LLt[, k] <- -loadings_LLt[, k]
    V_right[, k]      <- -V_right[, k]
  }
}

# CI propagation via fixed mean rotation
pca_lo <- pmin(L_q5 %*% V_right, L_q95 %*% V_right)
pca_hi <- pmax(L_q5 %*% V_right, L_q95 %*% V_right)

# ── PCA of Sigma_a = LL' + diag(psi^2) (total player) ─────────────────────────

Sigma_a       <- tcrossprod(L_mean) + diag(psi_a^2)
eg_Sa         <- eigen(Sigma_a, symmetric = TRUE)
var_Sa_full   <- eg_Sa$values                   # length P
var_share_Sa  <- var_Sa_full[1:K] / sum(var_Sa_full)  # first-K share of TOTAL Sigma_a trace

cat(sprintf("\nLambdaLambda' variance shares (K=%d):\n", K))
for (k in 1:K)
  cat(sprintf("  PC%d: %.1f%%  (eigenvalue %.4f)\n",
              k, 100 * var_share_LLt[k], var_LLt[k]))
cat(sprintf("\nSigma_a variance shares (first %d of %d PCs):\n", K, P))
for (k in 1:K)
  cat(sprintf("  PC%d: %.1f%%\n", k, 100 * var_share_Sa[k]))

# ── Write variance shares (append across K) ───────────────────────────────────

groups <- as.data.frame(feature_group_lookup())

var_tbl <- tibble(
  K            = K,
  pc           = paste0("PC", 1:K),
  eigenval_LLt = var_LLt,
  share_LLt    = var_share_LLt,
  eigenval_Sa  = var_Sa_full[1:K],
  share_Sa     = var_share_Sa
)
write_csv(var_tbl, file.path(paths$tables, sprintf("31b_k%d_variance_shares.csv", K)))
message(sprintf("Saved: 31b_k%d_variance_shares.csv", K))

# ── Build CI-table and reliable-count ─────────────────────────────────────────

pca_ci_tbl <- do.call(rbind, lapply(1:K, function(k) {
  tibble(
    K       = K,
    pc      = paste0("PC", k),
    feature = feature_names,
    group   = groups$group[match(feature_names, groups$feature)],
    loading = loadings_LLt[, k],
    ci_lo   = pca_lo[, k],
    ci_hi   = pca_hi[, k]
  ) |>
    mutate(reliable = (ci_lo > 0 & ci_hi > 0) | (ci_lo < 0 & ci_hi < 0))
}))
write_csv(pca_ci_tbl,
          file.path(paths$tables, sprintf("31b_k%d_pca_loading_ci.csv", K)))

rel_tbl <- pca_ci_tbl |>
  group_by(K, pc) |>
  summarise(
    n_features     = n(),
    n_reliable     = sum(reliable, na.rm = TRUE),
    share_reliable = mean(reliable, na.rm = TRUE),
    max_abs_load   = max(abs(loading), na.rm = TRUE),
    .groups        = "drop"
  )
write_csv(rel_tbl,
          file.path(paths$tables, sprintf("31b_k%d_reliable_counts.csv", K)))
message(sprintf("Saved: 31b_k%d_reliable_counts.csv", K))

cat("\nReliable-loader summary:\n")
print(rel_tbl)

# ── Append to cross-K comparison tables (idempotent) ──────────────────────────

var_all_path <- file.path(paths$tables, "31b_kcompare_variance_shares.csv")
rel_all_path <- file.path(paths$tables, "31b_kcompare_reliable_counts.csv")

update_comparison <- function(path, new_row) {
  existing <- if (file.exists(path)) {
    read_csv(path, show_col_types = FALSE) |> filter(K != !!K)
  } else NULL
  bind_rows(existing, new_row) |> arrange(K, pc) |> write_csv(path)
}
update_comparison(var_all_path, var_tbl)
update_comparison(rel_all_path, rel_tbl)
message("Appended into 31b_kcompare_variance_shares.csv and 31b_kcompare_reliable_counts.csv")

message(sprintf("\nStage 31b (K=%d) complete.", K))
