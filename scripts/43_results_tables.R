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

diag_raw <- read_csv(file.path(paths$tables, "28_diagnostics_cache.csv"), show_col_types = FALSE) |>
  mutate(family = sub("\\..*$", "", param))

# Exclude the LT-PD structural zeros (fixed-at-0 upper-triangular Lambda_a entries,
# 3 of them at K=3: [1,2],[1,3],[2,3]) -- zero variance across draws makes Rhat
# undefined (NaN) for these, not a real mixing problem. Also pull LLt.1.10 out as
# its own bellwether row rather than an n=1 "family".
diag <- diag_raw |> filter(family != "LLt", !is.na(rhat))
llt_diag <- diag_raw |> filter(family == "LLt")

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

diag_display <- diag_summary |>
  transmute(
    Parameter = family, N = n,
    `Median ESS(bulk)` = round(median_ess_bulk, 0),
    `Min ESS(bulk)` = round(min_ess_bulk, 0),
    `Max Rhat` = round(max_rhat, 3),
    `N(Rhat>1.01)` = n_rhat_gt_1.01
  ) |>
  bind_rows(tibble(
    Parameter = "LLt[1,10] (bellwether)", N = 1,
    `Median ESS(bulk)` = round(llt_diag$ess_bulk, 0), `Min ESS(bulk)` = round(llt_diag$ess_bulk, 0),
    `Max Rhat` = round(llt_diag$rhat, 3), `N(Rhat>1.01)` = as.integer(llt_diag$rhat > 1.01)
  ))
write_latex_table(diag_display, "43_diagnostics_summary")
message(sprintf("LLt[1,10] mixing bellwether (rotation-invariant, the goals x passes pair flagged in earlier sessions): ESS_bulk=%.0f, Rhat=%.4f -> %s",
                llt_diag$ess_bulk, llt_diag$rhat, if (llt_diag$ess_bulk > 200 && llt_diag$rhat < 1.01) "PASS" else "borderline"))

# ── §5.2 Rank-selection Criteria A/B/C recap (K=2,3,4) ────────────────────────

kc_resid <- read_csv(file.path(paths$tables, "29b_kcompare_summary.csv"), show_col_types = FALSE)
kc_var   <- read_csv(file.path(paths$tables, "31b_kcompare_variance_shares.csv"), show_col_types = FALSE)
kc_rel   <- read_csv(file.path(paths$tables, "31b_kcompare_reliable_counts.csv"), show_col_types = FALSE)

last_pc_per_k <- function(df, k) df |> filter(K == k) |> slice_max(pc, n = 1, with_ties = FALSE)

rank_recap <- lapply(c(2, 3, 4), function(k) {
  v <- last_pc_per_k(kc_var, k)
  rl <- last_pc_per_k(kc_rel, k)
  r  <- kc_resid |> filter(K == k)
  tibble(
    K = k,
    `Last-PC share of LLt` = sprintf("%.1f%%", 100 * v$share_LLt),
    `Last-PC reliable loaders` = sprintf("%d/%d", rl$n_reliable, rl$n_features),
    `Mean resid. > tau99` = sprintf("%.2f%%", 100 * r$mean_frac_gt_tau99),
    `nu_hat` = round(r$nu_hat, 2)
  )
})
rank_recap_tbl <- bind_rows(rank_recap)
write_csv(rank_recap_tbl, file.path(paths$tables, "43_rank_criteria_abc.csv"))
write_latex_table(rank_recap_tbl, "43_rank_criteria_abc", digits = 2)

# ── §5.3 Minutes-scaling (phi) diagnostic ─────────────────────────────────────

get_row <- function(path, var) {
  read_csv(path, show_col_types = FALSE) |> filter(variable == var)
}
nu_28  <- get_row(file.path(paths$tables, "28_real_lowrank_a_diag_b_t_k3_posterior_summary.csv"), "nu")
nu_34a <- get_row(file.path(paths$tables, "34a_real_lowrank_a_diag_b_t_mv_k3_posterior_summary.csv"), "nu")
nu_34b <- get_row(file.path(paths$tables, "34b_real_lowrank_a_diag_b_t_mv_phi_k3_posterior_summary.csv"), "nu")
phi_34b <- get_row(file.path(paths$tables, "34b_real_lowrank_a_diag_b_t_mv_phi_k3_posterior_summary.csv"), "phi")

phi_tbl <- tibble(
  Model = c("Constant sigma_e (stage 28)", "Minutes-scaled, phi=0.5 fixed (34a)", "Minutes-scaled, phi estimated (34b)"),
  `nu (90% CI)` = c(
    sprintf("%.2f [%.2f, %.2f]", nu_28$mean, nu_28$q5, nu_28$q95),
    sprintf("%.2f [%.2f, %.2f]", nu_34a$mean, nu_34a$q5, nu_34a$q95),
    sprintf("%.2f [%.2f, %.2f]", nu_34b$mean, nu_34b$q5, nu_34b$q95)
  ),
  `phi (90% CI)` = c("fixed at 0", "fixed at 0.5", sprintf("%.3f [%.3f, %.3f]", phi_34b$mean, phi_34b$q5, phi_34b$q95)),
  `Low-minute worst-cell share` = c("48.9% (88/180)", "8.3% (14/168)", "8.3% (14/168)")
)
write_csv(phi_tbl, file.path(paths$tables, "43_phi_comparison.csv"))
write_latex_table(phi_tbl, "43_phi_comparison")

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

message("Stage 43 complete.")
