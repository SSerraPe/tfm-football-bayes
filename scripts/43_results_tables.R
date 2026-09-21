# Stage 43 — Remaining native LaTeX table fragments for the Results chapter.
#
# Everything here reads already-existing small CSVs (no new CmdStan extraction) --
# unlike stages 40-42, this is pure reshaping + kableExtra formatting.
#
# Produces:
#   §5.1 model diagnostics       -> 43_diagnostics_summary.{csv,tex}
#   §5.2 rank-selection recap    -> 43_rank_criteria_abc.{csv,tex}
#   §5.3 minutes-scaling (phi)   -> 43_phi_comparison.{csv,tex}
#   §5.5 canonical loadings      -> 43_canonical_loadings_top.tex (per PC top loaders)
#   §5.6 archetype top/bottom    -> 43_archetype_pc{1,2,3}.tex

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/43_results_tables.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
source(file.path(model_root, "src", "latex_table_helpers.R"))
check_packages(c(required_base_packages))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

# ── §5.1 Model diagnostics ────────────────────────────────────────────────────

diag_raw <- read_csv(file.path(paths$tables, "34a_k3_diagnostics_cache.csv"), show_col_types = FALSE) |>
  mutate(family = sub("\\..*$", "", param))

# Item 1.1: "LLt" (all P*(P+1)/2 unique entries of Lambda_a Lambda_a') and "LLtEig" (the
# K eigenvalues d_1..d_K of the same matrix) are now full families, produced by
# B_block_b_diagnostics_28.R, and are summarised in exactly the same format as every
# other parameter family below -- no more single-cell bellwether.
# Exclude the LT-PD structural zeros (fixed-at-0 upper-triangular Lambda_a entries,
# 3 of them at K=3: [1,2],[1,3],[2,3]) -- zero variance across draws makes Rhat
# undefined (NaN) for these, not a real mixing problem.
diag <- diag_raw |> filter(!is.na(rhat))

diag_summary <- diag |>
  group_by(family) |>
  summarise(
    n = n(),
    median_ess_bulk = median(ess_bulk),
    min_ess_bulk = min(ess_bulk),
    max_rhat = max(rhat),
    n_rhat_gt_1.01 = sum(rhat > 1.01),
    .groups = "drop"
  ) |>
  arrange(family)
write_csv(diag_summary, file.path(paths$tables, "43_diagnostics_summary.csv"))

param_symbol <- c(
  Lambda_a = "$\\boldsymbol\\Lambda_a$", nu = "$\\nu$", psi_a = "$\\boldsymbol\\psi_a$",
  sigma_b = "$\\boldsymbol\\sigma_b$", sigma_e = "$\\boldsymbol\\sigma_e$",
  LLt = "$\\boldsymbol\\Lambda_a\\boldsymbol\\Lambda_a^\\top$ (all entries)",
  LLtEig = "$d_1,d_2,d_3$ (eigenvalues)"
)

diag_display <- diag_summary |>
  transmute(
    Parameter = param_symbol[family], `$N$` = n,
    `ESS$_\\text{bulk}$: median (min)` = sprintf("%.0f (%.0f)", median_ess_bulk, min_ess_bulk),
    `Max $\\hat R$` = round(max_rhat, 3),
    `$N(\\hat R{>}1.01)$` = n_rhat_gt_1.01
  )
write_latex_table(diag_display, "43_diagnostics_summary", escape = FALSE)

llt_summary <- diag_summary |> filter(family == "LLt")
eig_summary <- diag_summary |> filter(family == "LLtEig")
message(sprintf("Lambda_a Lambda_a' full-matrix summary (%d entries): median(min) ESS_bulk=%.0f(%.0f), max Rhat=%.4f, N(Rhat>1.01)=%d",
                llt_summary$n, llt_summary$median_ess_bulk, llt_summary$min_ess_bulk,
                llt_summary$max_rhat, llt_summary$n_rhat_gt_1.01))
message(sprintf("Eigenvalues d_1..d_K summary: median(min) ESS_bulk=%.0f(%.0f), max Rhat=%.4f, N(Rhat>1.01)=%d",
                eig_summary$median_ess_bulk, eig_summary$min_ess_bulk, eig_summary$max_rhat, eig_summary$n_rhat_gt_1.01))

# ── §5.2 Rank-selection Criteria A/B recap (K=2,3,4) ──────────────────────────
# NOTE (2026-09, thesis revision Phase 3): the "reliable loaders" column (the
# criterion formerly labelled B, now dropped from the thesis argument per
# methodology.tex Sec. 8 -- it gave no reason to prefer K=3 over K=4) is no
# longer included in the table. The residual-adequacy criterion is relabelled
# Criterion B throughout, matching methodology.tex. Real math-mode column
# headers replace the earlier raw-ASCII ones (K, "Last-PC share of LLt",
# "Mean resid. > tau99", "nu_hat").
# Item 2.3 (revision pass 3): nu_hat column dropped -- a column whose own caption
# had to disclaim its relevance shouldn't be a column. nu_hat still matters (Criterion
# B's tau99 threshold is each fit's own t(nu_hat)-implied value), but that's a half
# sentence in the prose, not implied by a fourth column here. nu_hat per fit is still
# in 29b_kcompare_summary.csv for anyone who wants it.

kc_resid <- read_csv(file.path(paths$tables, "29b_kcompare_summary.csv"), show_col_types = FALSE)
kc_var   <- read_csv(file.path(paths$tables, "31b_kcompare_variance_shares.csv"), show_col_types = FALSE)

last_pc_per_k <- function(df, k) df |> filter(K == k) |> slice_max(pc, n = 1, with_ties = FALSE)

rank_recap <- lapply(c(2, 3, 4), function(k) {
  v <- last_pc_per_k(kc_var, k)
  r  <- kc_resid |> filter(K == k)
  tibble(
    `$K$` = k,
    `Last-PC share of $\\bar{\\boldsymbol\\Lambda}_a\\bar{\\boldsymbol\\Lambda}_a^\\top$` = sprintf("%.1f\\%%", 100 * v$share_LLt),
    `Mean $\\Pr(|z|>\\tau_{99})$` = sprintf("%.2f\\%%", 100 * r$mean_frac_gt_tau99)
  )
})
rank_recap_tbl <- bind_rows(rank_recap)
write_csv(rank_recap_tbl, file.path(paths$tables, "43_rank_criteria_abc.csv"))
write_latex_table(rank_recap_tbl, "43_rank_criteria_abc", digits = 2, escape = FALSE)

# ── §5.3 Minutes-scaling (phi) diagnostic ─────────────────────────────────────

get_row <- function(path, var) {
  read_csv(path, show_col_types = FALSE) |> filter(variable == var)
}
nu_28  <- get_row(file.path(paths$tables, "28_real_lowrank_a_diag_b_t_k3_posterior_summary.csv"), "nu")
nu_34a <- get_row(file.path(paths$tables, "34a_real_lowrank_a_diag_b_t_mv_k3_posterior_summary.csv"), "nu")
nu_34b <- get_row(file.path(paths$tables, "34b_real_lowrank_a_diag_b_t_mv_phi_k3_posterior_summary.csv"), "nu")
phi_34b <- get_row(file.path(paths$tables, "34b_real_lowrank_a_diag_b_t_mv_phi_k3_posterior_summary.csv"), "phi")

phi_tbl <- tibble(
  # Item 3.2: description alone, no internal stage identifiers -- those stay in the repo.
  Model = c("Constant $\\sigma_e$", "$\\phi=0.5$ fixed (production)", "$\\phi$ estimated"),
  `$\\nu$ (90\\% CI)` = c(
    sprintf("%.2f [%.2f, %.2f]", nu_28$mean, nu_28$q5, nu_28$q95),
    sprintf("%.2f [%.2f, %.2f]", nu_34a$mean, nu_34a$q5, nu_34a$q95),
    sprintf("%.2f [%.2f, %.2f]", nu_34b$mean, nu_34b$q5, nu_34b$q95)
  ),
  `$\\phi$ (90\\% CI)` = c("fixed at 0", "fixed at 0.5", sprintf("%.3f [%.3f, %.3f]", phi_34b$mean, phi_34b$q5, phi_34b$q95)),
  `Low-minute share` = c("48.9\\% (88/180)", "8.3\\% (14/168)", "8.3\\% (14/168)")
)
write_csv(phi_tbl, file.path(paths$tables, "43_phi_comparison.csv"))
write_latex_table(phi_tbl, "43_phi_comparison", escape = FALSE)

# ── §5.5 Canonical loadings: top reliable loaders per PC ──────────────────────

pca_ci <- read_csv(file.path(paths$tables, "31_pca_loading_ci.csv"), show_col_types = FALSE)
for (pc_name in sort(unique(pca_ci$pc))) {
  top <- pca_ci |>
    filter(pc == pc_name, reliable) |>
    mutate(abs_loading = abs(loading)) |>
    arrange(desc(abs_loading)) |>
    slice_head(n = 8) |>
    transmute(
      Feature = gsub("_", " ", gsub("^per90_|^rate_", "", feature)),
      Group = group,
      Loading = sprintf("%.2f [%.2f, %.2f]", loading, ci_lo, ci_hi)
    )
  write_latex_table(top, paste0("43_canonical_loadings_", tolower(pc_name)))
}

# ── §5.6 Archetypes: top/bottom players per PC ────────────────────────────────

top_pc <- read_csv(file.path(paths$tables, "29_top_players_by_pc.csv"), show_col_types = FALSE)
for (pc_name in sort(unique(top_pc$pc))) {
  sub <- top_pc |> filter(pc == pc_name)
  disp <- bind_rows(
    sub |> filter(direction == "top")    |> arrange(rank_within) |> slice_head(n = 8),
    sub |> filter(direction == "bottom") |> arrange(rank_within) |> slice_head(n = 8)
  ) |>
    transmute(Player = player_name, Direction = direction,
              Score = round(.data[[paste0(tolower(pc_name), "_score")]], 2))
  write_latex_table(disp, paste0("43_archetype_", tolower(pc_name)))
}

# ── §5.6b Archetypes: combined PC1/PC2/PC3 table (thesis revision Phase 5) ────
# Replaces the three separate top/bottom-8 tables above with one compact
# top/bottom-4-per-axis table read side by side. Selection rule: within each
# axis/direction, the highest-magnitude players that are recognisable
# footballers, never reaching past roughly the top 15-20 by magnitude to find
# one -- except PC3-top, where the actual top-ranked players are all genuinely
# obscure (first recognisable name is rank 14) and the user chose the honest,
# unadjusted top 4 rather than reaching for recognisability (2026-09 decision).
#
# NOTE: pick by (pc, direction, player_name) against 29_top_players_by_pc.csv,
# not by player_name alone against 29b_player_archetypes.csv -- several common
# Spanish surnames (e.g. "J. Rodríguez") collide across unrelated players in
# the full 1,529-player table, and a name-only join silently pulled in the
# wrong player's score for one entry during Phase 5 (caught by inspection
# before it reached the thesis).
# SIGN CORRECTION (found during the "three replacement sections" revision, documented in
# full in scripts/48_style_space_scatter.R): 29_outlier_study.R's eigendecomposition of
# Lambda_bar %*% t(Lambda_bar) is never sign-anchored, unlike stage 31's
# (31_pca_with_ci.R), which anchors each PC so the feature with the largest |loading| is
# positive -- the convention Remark 4.24/rem:sign specifies and the one the new loading
# heatmap and style-space scatter (stages 47, 48) are built on. Checked directly: PC1
# happens to already agree with that convention; PC2 and PC3 are its exact reflection (an
# elite passer like Xavi scores pc2_score = -7.40 under stage 29's raw sign, when a
# passing-heavy player should score strongly positive on "passing volume/quality (+) vs.
# aerial/physical (-)"). 29_top_players_by_pc.csv inherits that same unanchored sign, so
# its "top"/"bottom" labels for PC2 and PC3 are reversed relative to the corrected
# convention used everywhere else in this chapter. Fixed here by flipping the retrieved
# score's sign for PC2/PC3 and swapping which (old) direction feeds the corrected top vs.
# bottom half of the table; PC1 needs neither, having already agreed.
# The flip_sign compensation this block used to need is gone: 29_outlier_study.R now
# anchors PC2/PC3's sign itself (largest-magnitude loading positive, Remark rem:sign),
# fixed 2026-09-20 at the source rather than worked around here. 29_top_players_by_pc.csv's
# "top" direction is therefore genuinely the positive pole for every PC, matching the
# loading heatmap and style-space scatter directly -- no relabelling needed.
pick_pc <- function(pc_name, dir_name, names) {
  sub <- top_pc |> filter(pc == pc_name, direction == dir_name)
  score_col <- paste0(tolower(pc_name), "_score")
  out <- sub |>
    filter(player_name %in% names) |>
    distinct(player_name, .keep_all = TRUE) |>
    mutate(player_name = factor(player_name, levels = names)) |>
    arrange(player_name)
  tibble(Player = as.character(out$player_name), Score = round(out[[score_col]], 2))
}

# Re-verified against the corrected production fit (Execution step F/G), not assumed
# stable -- the top-4/bottom-4 by magnitude changed on two of the six poles:
#  - PC1 top: Raul Albiol is no longer in the top ~15 by magnitude; replaced by Victor
#    Ruiz (#8 by magnitude, a recognisable Spain/Valencia/Napoli centre-back), keeping
#    Mascherano/Pique/Umtiti (still present, reordered).
#  - PC2 top: Kroos drops to #7 by magnitude; the new #2, L. Messi, is itself a top-4
#    name by magnitude and clearly recognisable, so no reach further down the ranking
#    was needed. Kept: Xavi, Iniesta, J. Rodriguez.
# PC1 bottom, PC2 bottom, and both PC3 poles are unchanged (same four names each,
# ranks 1-4 by magnitude, only minor internal reordering).
pc1_top    <- c("J. Mascherano", "Gerard Piqué", "Víctor Ruíz", "S. Umtiti")
pc1_bottom <- c("Neymar", "Vinícius Júnior", "K. Mbappé", "L. Messi")
pc2_positive <- c("Xavi", "L. Messi", "Andrés Iniesta", "J. Rodríguez")
pc2_negative <- c("C. Stuani", "A. Budimir", "Y. En-Nesyri", "S. Okazaki")
pc3_positive <- c("L. Messi", "Xavi", "M. Pjanić", "Cristiano Ronaldo")
pc3_negative <- c("Iván Alejo", "Juan Iglesias", "O. Vranješ", "Houboulang Mendes")

combined <- bind_cols(
  pick_pc("PC1", "top", pc1_top),
  pick_pc("PC2", "top", pc2_positive),
  pick_pc("PC3", "top", pc3_positive),
  .name_repair = "unique"
) |> bind_rows(
  bind_cols(
    pick_pc("PC1", "bottom", pc1_bottom),
    pick_pc("PC2", "bottom", pc2_negative),
    pick_pc("PC3", "bottom", pc3_negative),
    .name_repair = "unique"
  )
)
write_csv(combined, file.path(paths$tables, "43_archetype_combined.csv"))
# WARNING: the tracked tfm_mesio_latex_2026/tables/43_archetype_combined.tex is hand-edited
# on top of this script's raw output -- it carries a multicolumn PC-name header, a
# hspace-adjusted column spec, and a Top/Bottom \midrule + \multicolumn divider that
# write_latex_table()'s plain kable() call does not produce (results.tex just \input{}s the
# fragment directly; it does not add this formatting itself, despite what an earlier version
# of this message claimed). Discovered 2026-09-20 when a helper path fix (latex_table_helpers.R)
# caused this call to actually reach tfm_mesio_latex_2026/ for the first time in a while and
# clobbered the hand formatting; player data was unaffected (identical), only the header/divider
# were lost, and were restored from git/backup. If the underlying player rankings ever change,
# regenerate via this script, then manually reapply the multicolumn header and Top/Bottom divider
# -- don't just use the raw output.
write_latex_table(combined,
  col_names = c("Player", "Score", "Player", "Score", "Player", "Score"),
  name = "43_archetype_combined")
message("Combined archetype table written (RAW kable output -- reapply the hand-formatted ",
        "multicolumn header and Top/Bottom divider before using; see comment above).")

message("Stage 43 complete.")
