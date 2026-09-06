# Quick sanity report on the K=4 kcompare row.
# Prints all three kcompare CSVs (should now have K=2, K=3, K=4 rows)
# and flags whether the K=4 numbers meet the expected "noise floor" prediction.

suppressPackageStartupMessages({
  library(readr); library(dplyr)
})
setwd("/Users/sserra/Documents/01_MESIO_UPC/TFM/model")

cat("=== 29b_kcompare_summary.csv ===\n")
resid <- read_csv("outputs/tables/29b_kcompare_summary.csv", show_col_types = FALSE)
print(as.data.frame(resid))
cat("\n")

cat("=== 31b_kcompare_variance_shares.csv ===\n")
vs <- read_csv("outputs/tables/31b_kcompare_variance_shares.csv", show_col_types = FALSE)
print(as.data.frame(vs))
cat("\n")

cat("=== 31b_kcompare_reliable_counts.csv ===\n")
rc <- read_csv("outputs/tables/31b_kcompare_reliable_counts.csv", show_col_types = FALSE)
print(as.data.frame(rc))
cat("\n")

# ── Verdicts ─────────────────────────────────────────────────────────────────

if (!4 %in% resid$K) {
  cat("!!! K=4 row NOT in 29b_kcompare_summary — 29b hasn't run at K=4 yet.\n")
} else {
  k4_resid <- resid[resid$K == 4, ]
  cat(sprintf("K=4 residual adequacy: mean_frac_gt_tau99 = %.2f%%  (vs 1%% floor)\n",
              100 * k4_resid$mean_frac_gt_tau99))
  cat(sprintf("K=4 nu_hat = %.3f, tau_99 = %.3f\n",
              k4_resid$nu_hat, k4_resid$tau_99))
}

if (!4 %in% vs$K) {
  cat("!!! K=4 row NOT in 31b_kcompare_variance_shares — 31b hasn't run at K=4 yet.\n")
} else {
  k4_vs <- vs[vs$K == 4, ]
  pc4 <- k4_vs[k4_vs$pc == "PC4", ]
  if (nrow(pc4) == 1L) {
    cat(sprintf("K=4 PC4: eigenval = %.3f, share of LL' = %.1f%%, share of Sigma_a = %.1f%%\n",
                pc4$eigenval_LLt, 100 * pc4$share_LLt, 100 * pc4$share_Sa))
    # Verdict against expected "noise floor" (loose thresh)
    if (100 * pc4$share_LLt < 5) {
      cat("  ==> PC4 at NOISE FLOOR. K=3 vs K=4 verdict: K=3 wins on parsimony.\n")
    } else if (100 * pc4$share_LLt < 10) {
      cat("  ==> PC4 is small but non-trivial. K=3 vs K=4 verdict: K=3 preferred but borderline.\n")
    } else {
      cat("  ==> PC4 has meaningful variance share. K=3 vs K=4 verdict: AMBIGUOUS — supervisor discussion warranted.\n")
    }
    # Compare to K=3's PC3
    k3_pc3 <- vs[vs$K == 3 & vs$pc == "PC3", ]
    if (nrow(k3_pc3) == 1L) {
      cat(sprintf("  Reference: at K=3, PC3 share of LL' = %.1f%%\n",
                  100 * k3_pc3$share_LLt))
    }
  }
}

if (!4 %in% rc$K) {
  cat("!!! K=4 row NOT in 31b_kcompare_reliable_counts — 31b hasn't run at K=4 yet.\n")
} else {
  k4_rc <- rc[rc$K == 4, ]
  pc4_rc <- k4_rc[k4_rc$pc == "PC4", ]
  if (nrow(pc4_rc) == 1L) {
    cat(sprintf("K=4 PC4 reliable-loader count: %d/%d (%.1f%%)\n",
                pc4_rc$n_reliable, pc4_rc$n_features, 100 * pc4_rc$share_reliable))
    if (pc4_rc$n_reliable < 15) {
      cat("  ==> PC4 is noise-dominated. K=3 vs K=4 verdict: K=3 wins.\n")
    } else if (pc4_rc$n_reliable < 30) {
      cat("  ==> PC4 has partial identification. Borderline.\n")
    } else {
      cat("  ==> PC4 has strong identification. AMBIGUOUS.\n")
    }
  }
}

cat("\n=== Done. ===\n")
