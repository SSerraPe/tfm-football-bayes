# Block B: Stage 28 comprehensive diagnostics
#
# Extracts Lambda_a (144 = 48×3), psi_a (48), sigma_b (48), sigma_e (48), nu (1)
# using Python column-pass; computes ESS_bulk and Rhat via posterior package.
# Reports: nu median/CI/ESS/Rhat; Lambda_a ESS/Rhat (all, especially LLt[1,10]);
#          psi_a/sigma_b/sigma_e summaries; Pareto-k from existing posterior_summary;
#          count of Rhat > 1.01.
# Caches to outputs/tables/34a_k3_diagnostics_cache.csv.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/B_block_b_diagnostics_28.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "posterior"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(posterior) })

py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")
# Was hardcoded to a prior session's scratchpad dir (no longer exists on disk); switched
# to tempdir() to match the pattern used by the other extraction scripts (e.g. stage 40).
SCRATCHPAD <- file.path(tempdir(), "block_b_diagnostics_extract")
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

# Repointed to the corrected production fit (Execution step F): stage 34a, minutes-scaled,
# fixed phi=0.5, K=3 -- Execution step D confirmed K*=3 unchanged under this family.
MODEL_ID <- "34a_real_lowrank_a_diag_b_t_mv_k3"
csv_files <- discover_cmdstan_csv_files(MODEL_ID)
if (length(csv_files) == 0) stop("No CSV files found for ", MODEL_ID)
message("Found ", length(csv_files), " chain CSVs for ", MODEL_ID)

# ── Column index discovery (one parse of chain 1 header) ────────────────────

find_col_indices <- function(csv_path, param_names) {
  con <- file(csv_path, "r")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0L) return(setNames(rep(NA_integer_, length(param_names)), param_names))
    if (!startsWith(line, "#")) {
      hdr <- strsplit(line, ",")[[1]]
      idx <- match(param_names, hdr)
      return(setNames(idx, param_names))
    }
  }
}

# Discover all Lambda_a, psi_a, sigma_b, sigma_e, nu columns
message("Parsing CSV header to find column positions...")
K <- 3L; P <- 48L

lambda_names <- paste0("Lambda_a.", rep(1:P, K), ".", rep(1:K, each = P))
psi_names    <- paste0("psi_a.", 1:P)
sigma_b_names <- paste0("sigma_b.", 1:P)
sigma_e_names <- paste0("sigma_e.", 1:P)
all_params   <- c("nu", lambda_names, psi_names, sigma_b_names, sigma_e_names)

col_idx <- find_col_indices(csv_files[1], all_params)
missing  <- names(col_idx)[is.na(col_idx)]
if (length(missing) > 0) message("WARNING: not found in header: ", paste(missing, collapse=", "))

# Keep only params that were found and their 0-based col indices
found_params <- col_idx[!is.na(col_idx)]
message("Found ", length(found_params), " of ", length(all_params), " parameters")

# ── Batch extraction: one Python call per chain ──────────────────────────────

extract_batch <- function(csv_path, chain_id, params_0idx, out_path) {
  args_pairs <- paste0(names(params_0idx), ":", params_0idx - 1L)
  cmd_args <- c(shQuote(py_script), shQuote(csv_path),
                as.character(chain_id), shQuote(out_path),
                args_pairs)
  ret <- system2("python3", args = cmd_args, stdout = FALSE, stderr = FALSE)
  if (ret != 0) stop("Python extraction failed for chain ", chain_id)
  read.csv(out_path, check.names = FALSE)
}

message("Extracting ", length(found_params), " params from ", length(csv_files), " chains...")
chain_dfs <- lapply(seq_along(csv_files), function(i) {
  out <- file.path(SCRATCHPAD, sprintf("b_block_b_chain%d.csv", i))
  message("  Chain ", i, "...")
  extract_batch(csv_files[i], i, found_params, out)
})

draws_df <- do.call(rbind, chain_dfs)
n_iter   <- max(draws_df[[".iteration"]])
n_chains <- max(draws_df[[".chain"]])
message(sprintf("Combined: %d iterations × %d chains = %d total draws",
                n_iter, n_chains, nrow(draws_df)))

# ── Build draws array and compute ESS/Rhat for each param ───────────────────

param_cols <- setdiff(names(draws_df), c(".chain", ".iteration"))

make_array <- function(df, param, n_iter, n_chains) {
  a <- array(NA_real_, dim = c(n_iter, n_chains, 1L))
  for (ch in seq_len(n_chains)) a[, ch, 1L] <- df[[param]][df[[".chain"]] == ch]
  as_draws_array(a)
}

message("Computing ESS and Rhat for ", length(param_cols), " parameters...")

diagnostics <- lapply(param_cols, function(p) {
  arr <- make_array(draws_df, p, n_iter, n_chains)
  data.frame(
    param     = p,
    mean      = mean(draws_df[[p]]),
    sd        = sd(draws_df[[p]]),
    q05       = quantile(draws_df[[p]], 0.05),
    q50       = quantile(draws_df[[p]], 0.50),
    q95       = quantile(draws_df[[p]], 0.95),
    ess_bulk  = ess_bulk(arr),
    rhat      = rhat(arr),
    row.names = NULL
  )
})
diag_df <- do.call(rbind, diagnostics)

# ── nu summary ───────────────────────────────────────────────────────────────

nu_row <- diag_df[diag_df$param == "nu", ]
message("\n=== nu ===")
message(sprintf("  median = %.3f  [%.3f, %.3f]  ESS = %.0f  Rhat = %.3f",
                nu_row$q50, nu_row$q05, nu_row$q95,
                nu_row$ess_bulk, nu_row$rhat))

# ── Lambda_a summary ─────────────────────────────────────────────────────────

la_diag <- diag_df[startsWith(diag_df$param, "Lambda_a."), ]

# Lower-tri mask: Lambda_a.i.k where k <= i
la_diag <- la_diag %>%
  mutate(
    row_i = as.integer(sub("Lambda_a\\.(\\d+)\\.(\\d+)", "\\1", param)),
    col_k = as.integer(sub("Lambda_a\\.(\\d+)\\.(\\d+)", "\\2", param)),
    is_lower_tri = (col_k <= row_i)
  )

la_lt <- la_diag[la_diag$is_lower_tri, ]  # 141 entries

message("\n=== Lambda_a (lower-triangular, ", nrow(la_lt), " entries) ===")
message(sprintf("  ESS_bulk: min = %.0f  median = %.0f  max = %.0f",
                min(la_lt$ess_bulk), median(la_lt$ess_bulk), max(la_lt$ess_bulk)))
message(sprintf("  Rhat:     min = %.3f  median = %.3f  max = %.3f",
                min(la_lt$rhat), median(la_lt$rhat), max(la_lt$rhat)))
worst_ess <- la_lt[which.min(la_lt$ess_bulk), ]
worst_rht <- la_lt[which.max(la_lt$rhat), ]
message(sprintf("  Worst ESS: %s (ESS=%.0f, Rhat=%.3f)",
                worst_ess$param, worst_ess$ess_bulk, worst_ess$rhat))
message(sprintf("  Worst Rhat: %s (ESS=%.0f, Rhat=%.3f)",
                worst_rht$param, worst_rht$ess_bulk, worst_rht$rhat))

# ── Item 1.1: full Lambda_a Lambda_a' convergence summary ────────────────────
# Replaces the single-cell LLt[1,10] bellwether with a summary over the FULL,
# rotation-invariant common-variance matrix C = Lambda_a %*% t(Lambda_a): all
# P*(P+1)/2 unique entries (P diagonal + P(P-1)/2 off-diagonal), computed draw-by-draw,
# plus the K eigenvalues d_1 >= ... >= d_K that Criterion A and the variance-share table
# actually depend on. Raw Lambda_a mixes less cleanly under the LT-PD identification
# constraint (Remark rem:ltpd_notinterp); this checks the object that is actually of
# substantive interest, in full, rather than by example.

lambda_cols        <- lambda_names
lambda_draws_by_row <- as.matrix(draws_df[, lambda_cols, drop = FALSE])
storage.mode(lambda_draws_by_row) <- "double"
n_draws_total <- nrow(lambda_draws_by_row)

pair_idx <- which(upper.tri(matrix(0, P, P), diag = TRUE), arr.ind = TRUE)
n_pairs  <- nrow(pair_idx)
stopifnot(n_pairs == P * (P + 1) / 2)

llt_mat <- matrix(NA_real_, n_draws_total, n_pairs)
eig_mat <- matrix(NA_real_, n_draws_total, K)

message("\nComputing full Lambda_a Lambda_a' (", n_pairs, " entries) and its top-", K,
        " eigenvalues for all ", n_draws_total, " draws...")
for (s in seq_len(n_draws_total)) {
  L_s <- matrix(lambda_draws_by_row[s, ], P, K)
  C_s <- tcrossprod(L_s)
  llt_mat[s, ] <- C_s[pair_idx]
  eig_mat[s, ] <- eigen(C_s, symmetric = TRUE, only.values = TRUE)$values[seq_len(K)]
  if (s %% 1000 == 0) message("  draw ", s, "/", n_draws_total)
}

chain_iter_df <- data.frame(.chain = draws_df[[".chain"]], .iteration = draws_df[[".iteration"]])
ess_rhat_of <- function(v) {
  arr <- make_array(cbind(chain_iter_df, v = v), "v", n_iter, n_chains)
  c(ess_bulk = ess_bulk(arr), rhat = rhat(arr))
}

llt_diag_stats <- t(vapply(seq_len(n_pairs), function(j) ess_rhat_of(llt_mat[, j]), numeric(2)))
eig_diag_stats <- t(vapply(seq_len(K),        function(k) ess_rhat_of(eig_mat[, k]), numeric(2)))

message(sprintf("\n=== Lambda_a Lambda_a' (full, %d unique entries) ===", n_pairs))
message(sprintf("  ESS_bulk: median = %.0f  min = %.0f", median(llt_diag_stats[, "ess_bulk"]), min(llt_diag_stats[, "ess_bulk"])))
message(sprintf("  Rhat:     max = %.4f", max(llt_diag_stats[, "rhat"])))
message(sprintf("  N(Rhat > 1.01): %d of %d (%.1f%%) -- with %d entries, a handful above 1.01 is expected by chance",
                sum(llt_diag_stats[, "rhat"] > 1.01), n_pairs, 100 * mean(llt_diag_stats[, "rhat"] > 1.01), n_pairs))

message(sprintf("\n=== Eigenvalues d_1..d_%d of Lambda_a Lambda_a' ===", K))
for (k in seq_len(K)) {
  message(sprintf("  d_%d: ESS=%.0f  Rhat=%.4f", k, eig_diag_stats[k, "ess_bulk"], eig_diag_stats[k, "rhat"]))
}

# Caution (b): Rhat is unreliable for entries concentrated near zero (many off-diagonal
# rank-3 entries near the LT-PD boundary). Report a magnitude-conditioned summary
# alongside the unconditional one so a genuinely misleading case isn't masked -- without
# dropping inconvenient cells from the unconditional summary above.
llt_post_mean <- colMeans(llt_mat)
above_floor   <- abs(llt_post_mean) > 0.05
message(sprintf("\n  Conditional on |posterior mean| > 0.05 (%d of %d entries): max Rhat = %.4f, N(Rhat>1.01) = %d",
                sum(above_floor), n_pairs, max(llt_diag_stats[above_floor, "rhat"]),
                sum(llt_diag_stats[above_floor, "rhat"] > 1.01)))

# Per-factor Lambda_a ESS
for (k in 1:K) {
  la_k <- la_lt[la_lt$col_k == k, ]
  message(sprintf("  Factor %d: ESS min=%.0f  median=%.0f  max=%.0f  | Rhat max=%.3f",
                  k, min(la_k$ess_bulk), median(la_k$ess_bulk), max(la_k$ess_bulk), max(la_k$rhat)))
}

# ── psi_a, sigma_b, sigma_e ──────────────────────────────────────────────────

for (grp in c("psi_a", "sigma_b", "sigma_e")) {
  g_df <- diag_df[startsWith(diag_df$param, paste0(grp, ".")), ]
  message(sprintf("\n=== %s (%d entries) ===", grp, nrow(g_df)))
  message(sprintf("  ESS_bulk: min=%.0f  median=%.0f  max=%.0f",
                  min(g_df$ess_bulk), median(g_df$ess_bulk), max(g_df$ess_bulk)))
  message(sprintf("  Rhat:     min=%.3f  median=%.3f  max=%.3f",
                  min(g_df$rhat), median(g_df$rhat), max(g_df$rhat)))
}

# ── Rhat > 1.01 count ────────────────────────────────────────────────────────

n_rhat_101 <- sum(diag_df$rhat > 1.01, na.rm = TRUE)
n_rhat_105 <- sum(diag_df$rhat > 1.05, na.rm = TRUE)
message(sprintf("\n=== Global Rhat summary ==="))
message(sprintf("  Rhat > 1.01: %d of %d parameters (%.1f%%)",
                n_rhat_101, nrow(diag_df), 100 * n_rhat_101 / nrow(diag_df)))
message(sprintf("  Rhat > 1.05: %d of %d parameters (%.1f%%)",
                n_rhat_105, nrow(diag_df), 100 * n_rhat_105 / nrow(diag_df)))

# ── Pareto-k from existing posterior summary ──────────────────────────────────

ps_path <- file.path(paths$tables, "28_real_lowrank_a_diag_b_t_k3_posterior_summary.csv")
if (file.exists(ps_path)) {
  ps <- read_csv(ps_path, show_col_types = FALSE)
  if ("khat" %in% names(ps) || "pareto_k" %in% names(ps)) {
    k_col <- if ("pareto_k" %in% names(ps)) "pareto_k" else "khat"
    k_vals <- ps[[k_col]]
    message(sprintf("\n=== Pareto-k (from %s) ===", basename(ps_path)))
    message(sprintf("  Share > 0.5: %.1f%%   > 0.7: %.1f%%",
                    100 * mean(k_vals > 0.5, na.rm = TRUE),
                    100 * mean(k_vals > 0.7, na.rm = TRUE)))
  } else {
    message("\nNote: no Pareto-k in existing posterior summary (LOO not yet run for stage 28)")
  }
} else {
  message("\nNo stage 28 posterior summary found for Pareto-k (LOO not yet run)")
}

# ── Save cache ────────────────────────────────────────────────────────────────
# Full Lambda_a Lambda_a' summary: one row per unique entry, named "LLt.<p>.<q>" (p<=q)
# so that 43_results_tables.R's family-prefix grouping ("LLt") aggregates them the same
# way it aggregates Lambda_a/psi_a/etc. Eigenvalue rows are named "LLtEig.<k>" (own family).

llt_rows <- data.frame(
  param    = sprintf("LLt.%d.%d", pair_idx[, "row"], pair_idx[, "col"]),
  mean     = colMeans(llt_mat), sd = apply(llt_mat, 2, sd),
  q05      = apply(llt_mat, 2, quantile, 0.05), q50 = apply(llt_mat, 2, quantile, 0.50),
  q95      = apply(llt_mat, 2, quantile, 0.95),
  ess_bulk = llt_diag_stats[, "ess_bulk"], rhat = llt_diag_stats[, "rhat"],
  row.names = NULL
)
eig_rows <- data.frame(
  param    = sprintf("LLtEig.%d", seq_len(K)),
  mean     = colMeans(eig_mat), sd = apply(eig_mat, 2, sd),
  q05      = apply(eig_mat, 2, quantile, 0.05), q50 = apply(eig_mat, 2, quantile, 0.50),
  q95      = apply(eig_mat, 2, quantile, 0.95),
  ess_bulk = eig_diag_stats[, "ess_bulk"], rhat = eig_diag_stats[, "rhat"],
  row.names = NULL
)
cache <- rbind(diag_df, llt_rows, eig_rows)

write_csv(cache, file.path(paths$tables, "34a_k3_diagnostics_cache.csv"))
message("\n34a_k3_diagnostics_cache.csv saved (", nrow(cache), " rows)")
message("Block B complete.")
