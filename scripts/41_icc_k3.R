# Stage 41 — Per-feature ICC table for the production fit (stage 34a, minutes-scaled,
# fixed phi=0.5, K*=3 -- repointed from the constant-scale stage 28 in revision pass 3,
# Execution step F, after Execution step D confirmed K*=3 unchanged).
#
# Chapter 4 (methodology.tex, rem:bayesian_icc, def:icc) promises: "Posterior medians and
# 95% credible intervals are reported in Chapter~\ref{chap:results}" for ICC_player,
# ICC_season, ICC_residual per feature, evaluated at the reference exposure m_n=m_ref
# (Decision 3 of the revision that corrected this chapter). These are exactly the
# `prop_a`/`prop_b`/`prop_e` generated quantities already present in the posterior summary
# CSV (confirmed: no new CmdStan-CSV extraction needed) -- confirmed by direct reading of
# the Stan file that these are already computed at m_ref (var_e uses the base, unscaled
# sigma_e[p]), so no change to this extraction logic is needed for that decision, only the
# input file.
#
# Reads:  outputs/tables/34a_real_lowrank_a_diag_b_t_mv_k3_posterior_summary.csv
#         data/processed/model_objects.rds (feature name order)
# Writes: outputs/tables/41_icc_summary_k3.csv
#         tfm_mesio_latex_2026/tables/41_icc_top15.tex   (main-text table: top 15 by ICC_player)
#         tfm_mesio_latex_2026/tables/41_icc_full48.tex  (appendix-style full 48-row longtable body)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/41_icc_k3.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
source(file.path(model_root, "src", "latex_table_helpers.R"))
check_packages(c(required_base_packages))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

mo <- readRDS(file.path(paths$processed, "model_objects.rds"))
feature_names <- mo$variable_names
P <- length(feature_names)

ps_path <- file.path(paths$tables, "34a_real_lowrank_a_diag_b_t_mv_k3_posterior_summary.csv")
ps <- read_csv(ps_path, show_col_types = FALSE)

extract_prop <- function(prefix) {
  ps |>
    filter(grepl(paste0("^", prefix, "\\["), variable)) |>
    mutate(p = as.integer(sub(paste0(prefix, "\\[(\\d+)\\]"), "\\1", variable))) |>
    arrange(p) |>
    transmute(p, median = median, q5 = q5, q95 = q95)
}

prop_a <- extract_prop("prop_a")
prop_b <- extract_prop("prop_b")
prop_e <- extract_prop("prop_e")
stopifnot(nrow(prop_a) == P, nrow(prop_b) == P, nrow(prop_e) == P)

icc_tbl <- tibble(
  feature       = feature_names,
  icc_player     = prop_a$median, icc_player_q5  = prop_a$q5, icc_player_q95  = prop_a$q95,
  icc_season     = prop_b$median, icc_season_q5  = prop_b$q5, icc_season_q95  = prop_b$q95,
  icc_residual   = prop_e$median, icc_residual_q5 = prop_e$q5, icc_residual_q95 = prop_e$q95
) |>
  arrange(desc(icc_player))

write_csv(icc_tbl, file.path(paths$tables, "41_icc_summary_k3.csv"))
message("Saved: 41_icc_summary_k3.csv")

fmt_ci <- function(m, lo, hi) sprintf("%.2f [%.2f, %.2f]", m, lo, hi)

display_tbl <- icc_tbl |>
  transmute(
    Feature       = gsub("_", " ", gsub("^per90_|^rate_", "", feature)),
    `ICC(player)` = fmt_ci(icc_player, icc_player_q5, icc_player_q95),
    `ICC(season)` = fmt_ci(icc_season, icc_season_q5, icc_season_q95),
    `ICC(resid.)` = fmt_ci(icc_residual, icc_residual_q5, icc_residual_q95)
  )

# WARNING: 41_icc_full48.tex is hand-edited on top of this script's raw output -- it
# carries a \caption{}/\label{}/\endfirsthead/\endhead longtable header block
# (write_latex_table() does not support inline longtable captions) so it gets a real,
# auto-numbered table float. Back up the current file, run this script, then manually
# reapply that header/footer wrapper around the fresh data rows -- do not use the raw
# output directly. (This bit for real 2026-09-20 when latex_table_helpers.R's output path
# was fixed to actually reach tfm_mesio_latex_2026/ instead of a defunct tfm_latex/ copy;
# "41_icc_top15.tex" is unused in the main body, replaced by stage 45's stacked-bar figure.)
write_latex_table(head(display_tbl, 15), "41_icc_top15")
write_latex_table(display_tbl, "41_icc_full48", longtable = TRUE)

message("Stage 41 complete.")
