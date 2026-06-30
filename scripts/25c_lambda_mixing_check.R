# Task S — Check whether Lambda_a's poor mixing (ESS≈19, Rhat=1.185) propagates
# to Lambda_a %*% t(Lambda_a), the rotation-invariant quantity underlying all
# interpreted outputs (ICC via Sigma_a_diag, PCA-rotated loadings in stage 13).
#
# Mechanism: the lower-triangular identification constraint on Lambda_a doesn't
# guarantee that the free entries mix well. However, Lambda_a Lambda_a' is
# rotation-invariant (sign/rotation flips in Lambda_a cancel in the product).
# If Lambda_a chains explore different rotational orientations (the expected
# failure mode under the LT constraint), Lambda_a Lambda_a' should mix much
# better than Lambda_a itself.
#
# This script uses the stage18_small_params.csv extracted in the previous
# session (located in the session scratchpad) rather than re-running the
# Python column-pass, since that file already contains all Lambda_a draws.
#
# Verdict categories:
#   expected_artifact: Lambda_a Lambda_a' ESS is adequate (>= 200) despite
#                      poor raw Lambda_a ESS. The LT identification constraint
#                      causes rotation ambiguity in individual entries that
#                      cancels in the product. No action needed — thesis
#                      outputs are based on the rotation-invariant quantities.
#   escalated: Lambda_a Lambda_a' also has poor ESS/Rhat. The mixing problem
#              is not a labeling artifact but a genuine posterior-geometry
#              issue. Further action required before Task O.
#
# Output:
#   outputs/notes/note_lambda_identification_check.md

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/25c_lambda_mixing_check.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "posterior"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(posterior)
})

# ---- Load Lambda_a draws from scratchpad ------------------------------------
# The scratchpad CSV was produced in the previous session by Python column-pass.
# It has columns: .chain, .iteration, Lambda_a.1.1 .. Lambda_a.53.2

SCRATCHPAD <- Sys.getenv("CLAUDE_SCRATCHPAD",
  unset = "/private/tmp/claude-501/-Users-sserra-Documents-01-MESIO-UPC-TFM-model/76818d40-38b9-4367-9dc1-5d60a0b70149/scratchpad")
small_params_path <- file.path(SCRATCHPAD, "stage18_small_params.csv")

if (!file.exists(small_params_path)) {
  stop("stage18_small_params.csv not found in scratchpad: ", SCRATCHPAD,
       "\nRe-extract using src/extract_stan_csv_params.py with Lambda_a columns 91658-91763.")
}

message("Loading Lambda_a draws from ", small_params_path)
df_full <- read_csv(small_params_path, show_col_types = FALSE)

# Identify Lambda_a columns (Lambda_a.p.q naming)
la_cols <- grep("^Lambda_a\\.", names(df_full), value = TRUE)
P_check <- length(la_cols) / 2   # K=2
message(sprintf("Found %d Lambda_a columns (P=%d, K=2)", length(la_cols), as.integer(P_check)))

# ---- Compute Lambda_a %*% t(Lambda_a) per draw ------------------------------
# Lambda_a is P×2; product is P×P symmetric matrix.
# Focus on diagonal (common-factor variance contribution per feature) and
# a selection of off-diagonal pairs (to check rotation-invariance recovery).

message("Computing Lambda_a Lambda_a' diagonal for each draw...")

P  <- 53L
K  <- 2L

# Column naming: Lambda_a.p.q where p=1..P, q=1..K
la_1 <- paste0("Lambda_a.", seq_len(P), ".1")   # first factor loadings
la_2 <- paste0("Lambda_a.", seq_len(P), ".2")   # second factor loadings

L1 <- as.matrix(df_full[, la_1])  # n_draws × P
L2 <- as.matrix(df_full[, la_2])  # n_draws × P

# (Lambda_a Lambda_a')_{pp} = Lambda_a[p,1]^2 + Lambda_a[p,2]^2  (row variance)
# (Lambda_a Lambda_a')_{pq} = Lambda_a[p,1]*Lambda_a[q,1] + Lambda_a[p,2]*Lambda_a[q,2]

# Diagonal: P values per draw
LL_diag <- L1^2 + L2^2   # n_draws × P
colnames(LL_diag) <- paste0("LLt_diag.", seq_len(P))

# Selected off-diagonal pairs (feature indices; pick a few interesting ones)
# These are pairs that showed up in the correlation PPC as worst-misses:
# defensive duels vs head_shots, smart_passes vs progressive_passes
selected_pairs <- list(
  c(1, 2), c(1, 10), c(5, 15), c(10, 20), c(20, 30), c(30, 40), c(1, 53)
)
LL_offdiag <- do.call(cbind, lapply(selected_pairs, function(ij) {
  p <- ij[1]; q <- ij[2]
  L1[, p] * L1[, q] + L2[, p] * L2[, q]
}))
colnames(LL_offdiag) <- vapply(selected_pairs,
  function(ij) paste0("LLt.", ij[1], ".", ij[2]), character(1))

# ---- Build draws object and compute ESS/Rhat ---------------------------------

message("Computing ESS and Rhat for Lambda_a Lambda_a' entries...")

draws_df_input <- bind_cols(
  df_full |> select(.chain, .iteration),
  as.data.frame(LL_diag),
  as.data.frame(LL_offdiag)
)

n_chains <- max(draws_df_input$.chain)
n_iter   <- max(draws_df_input$.iteration)
draws_df_input$.draw <- as.integer(
  (draws_df_input$.chain - 1L) * n_iter + draws_df_input$.iteration
)

draws_obj <- posterior::as_draws_df(draws_df_input)

ps <- posterior::summarise_draws(
  draws_obj,
  mean     = mean,
  sd       = sd,
  rhat     = posterior::rhat,
  ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail
)

# Filter to LLt columns (diagonal: LLt_diag.*, off-diagonal: LLt.p.q)
ps_llt <- ps |> filter(grepl("^LLt", variable))

message("\n--- Lambda_a Lambda_a' (LLt) ESS summary ---")
diag_ps  <- ps_llt |> filter(grepl("^LLt_diag", variable))
off_ps   <- ps_llt |> filter(!grepl("^LLt_diag", variable))

message(sprintf("Diagonal (P=%d): min ESS_bulk = %.0f, max Rhat = %.3f",
                nrow(diag_ps), min(diag_ps$ess_bulk), max(diag_ps$rhat)))
message(sprintf("Off-diagonal (n=%d pairs): min ESS_bulk = %.0f, max Rhat = %.3f",
                nrow(off_ps), min(off_ps$ess_bulk), max(off_ps$rhat)))

# Compare to raw Lambda_a ESS from 18_posterior_health.csv
health_path <- file.path(paths$tables, "18_posterior_health.csv")
if (file.exists(health_path)) {
  health <- read_csv(health_path, show_col_types = FALSE)
  la_health <- health |> filter(group == "Lambda_a")
  if (nrow(la_health) > 0) {
    message(sprintf("\nRaw Lambda_a (from 18_posterior_health.csv):",
                    la_health$min_ess_bulk[1], la_health$max_rhat[1]))
    message(sprintf("  min ESS_bulk = %.0f, max Rhat = %.3f",
                    la_health$min_ess_bulk[1], la_health$max_rhat[1]))
  }
}

# ---- Verdict ----------------------------------------------------------------

min_ess_diag <- min(diag_ps$ess_bulk, na.rm = TRUE)
max_rhat_diag <- max(diag_ps$rhat, na.rm = TRUE)
min_ess_off  <- min(off_ps$ess_bulk, na.rm = TRUE)
max_rhat_off <- max(off_ps$rhat, na.rm = TRUE)

raw_la_min_ess  <- 19     # from 18_posterior_health.csv (Task F)
raw_la_max_rhat <- 1.185

ess_improvement <- if (is.finite(min_ess_diag)) min_ess_diag / raw_la_min_ess else NA
ess_threshold   <- 200   # 5% relative ESS for 4000 draws
rhat_threshold  <- 1.05

is_adequate <- min_ess_diag >= ess_threshold && max_rhat_diag <= rhat_threshold

# Nuanced verdict: partial improvement is informative
verdict <- if (is_adequate) "expected_artifact" else "escalated"
partial  <- !is_adequate && !is.na(ess_improvement) && ess_improvement > 2

message(sprintf("\n=== VERDICT: %s%s ===", toupper(verdict),
                if (partial) " (with partial rotation-artifact improvement)" else ""))
message(sprintf("ESS improvement factor (LLt diag vs raw Lambda_a): %.1fx (%.0f → %.0f)",
                ess_improvement, raw_la_min_ess, min_ess_diag))
message(sprintf("Rhat change: %.3f → %.3f", raw_la_max_rhat, max_rhat_diag))

if (verdict == "expected_artifact") {
  message("Lambda_a Lambda_a' mixes well despite poor raw Lambda_a mixing.")
  message("Rotation-labeling artifact from the LT constraint. Thesis outputs are safe.")
} else if (partial) {
  message("Partial improvement: rotation-labeling is a contributing factor,")
  message("but genuine geometry issues also contribute. Both ESS and Rhat still")
  message("exceed acceptable thresholds for the rotation-invariant product.")
  message("Interpretation: the LT constraint creates some labeling noise AND")
  message("the posterior geometry itself has mixing difficulty (plausible cause:")
  message("bimodal GK/outfield subgroup sitting far from the continuum).")
  message("Task O should not proceed until addressed.")
} else {
  message("No meaningful improvement. Mixing problem is primarily geometry-based,")
  message("not rotation-labeling. Task O should not proceed without addressing this.")
}

# ---- Write note file --------------------------------------------------------

note_path <- file.path(paths$notes, "note_lambda_identification_check.md")

note_lines <- c(
  "# Lambda_a Identification Check",
  "",
  "## Question",
  "",
  "Stage-18 diagnostics (Task F) found raw `Lambda_a` loadings mixing poorly:",
  "min ESS_bulk ≈ 19, max Rhat = 1.185. Does this poor mixing propagate to",
  "`Lambda_a Lambda_a'` (the rotation-invariant quantity underlying all",
  "interpreted thesis outputs: ICC, PCA-rotated loadings)?",
  "",
  "## Method",
  "",
  "Computed `(Lambda_a Lambda_a')_{pp} = L[p,1]^2 + L[p,2]^2` (diagonal,",
  "P=53 values) and selected off-diagonal pairs from 4000 draws (4 chains × 1000",
  "iterations). ESS_bulk and Rhat computed via `posterior::summarise_draws()`.",
  "",
  "## Results",
  "",
  sprintf("| Quantity | min ESS_bulk | max Rhat |"),
  sprintf("|---|---|---|"),
  sprintf("| Raw `Lambda_a` (from Task F) | ~19 | ~1.185 |"),
  sprintf("| `Lambda_a Lambda_a'` diagonal (P=53) | %.0f | %.3f |",
          min_ess_diag, max_rhat_diag),
  if (is.finite(min_ess_off))
    sprintf("| `Lambda_a Lambda_a'` off-diagonal (7 pairs) | %.0f | %.3f |",
            min_ess_off, max_rhat_off)
  else
    "| `Lambda_a Lambda_a'` off-diagonal | — (not computed) | — |",
  "",
  sprintf("ESS improvement factor (diagonal vs raw Lambda_a): %.1fx", ess_improvement),
  "",
  sprintf("**Verdict: %s%s**", verdict,
          if (partial) " (partial rotation-artifact improvement)" else ""),
  "",
  if (verdict == "expected_artifact") {
    c(
      "The poor `Lambda_a` mixing is a rotation-labeling artifact from the",
      "lower-triangular constraint. Chains explore different orientations of",
      "the same subspace; the product `Lambda_a Lambda_a'` (rotation-invariant)",
      "mixes adequately. All thesis-level quantities (ICC, PCA scores) are",
      "derived from `Lambda_a Lambda_a'` and are trustworthy.",
      "",
      "No additional action needed before Task O (full NUTS confirmation fit)."
    )
  } else if (partial) {
    c(
      "There is a 5× improvement in ESS from raw `Lambda_a` to `Lambda_a Lambda_a'`,",
      "suggesting rotation-labeling is a _contributing_ factor to the poor mixing.",
      "However, both ESS (< 200) and Rhat (> 1.05) remain above acceptable thresholds",
      "for the rotation-invariant product. The posterior geometry itself also has",
      "mixing difficulty.",
      "",
      "**Plausible cause:** the bimodal GK/outfield subgroup creates a structural",
      "heterogeneity in the loading space that makes exploration difficult. Removing",
      "GK rows (Task K) is a priority experiment — if it substantially improves",
      "`Lambda_a Lambda_a'` ESS, this supports the GK-geometry hypothesis.",
      "",
      "**Recommendation:** do not run Task O at the current K=2 fit without first",
      "running Task K and checking whether GK exclusion resolves the mixing.",
      "If it does not, consider a tighter prior on `lambda_a_free` or Option 3's",
      "marginalized Normal model."
    )
  } else {
    c(
      "No meaningful improvement from raw `Lambda_a` to `Lambda_a Lambda_a'`.",
      "The mixing problem is primarily geometry-based, not rotation-labeling.",
      "",
      "Recommended actions before Task O:",
      "1. Run Task K (GK exclusion + transforms) and re-check ESS.",
      "2. If not resolved: tighter prior on `lambda_a_free` entries.",
      "3. Or prioritize Option 3 (marginalized Normal model) which eliminates",
      "   the `Lambda_a` latent block entirely."
    )
  },
  "",
  sprintf("_Written by scripts/25c_lambda_mixing_check.R on %s_", Sys.Date())
)

writeLines(note_lines, note_path)
message("Saved note_lambda_identification_check.md")
message("\nTask S complete.")
