# SUPERSEDED / PROTOTYPE (migrated 2026-09-06): standalone port of the Bayesian
# Multidimensional Scaling (BMDS) exploration originally in
# tfm_playground/playground_BMDS.Rmd, plus the classical CMDS/PCA/FA sanity checks
# that preceded it in that notebook. This is NOT part of the active thesis pipeline:
# no current thesis chapter reports BMDS results (Chapter 3 covers only classical
# PCA/FA/CMDS; the production Bayesian model in Chapters 4-5 is the additive
# player+season+residual factor model, not BMDS). It is preserved here for provenance,
# since a completed BMDS run exists and could inform future extensions of this thesis.
#
# Dependencies: this script requires the `bayMDS`, `psych`, and `ggrepel` packages.
#
# The completed MCMC fit this script loads (`bmds_cosine_subset_N900_p2_*.rds`, ~2.9GB)
# is NOT copied into the tracked project structure -- it stays under
# tfm_playground/artifacts/, referenced below by relative path from the model/ root.
# Re-running bmdsMCMC() from scratch on a fresh subset is possible (the call is left
# in place below, guarded by `RUN_BMDS_MCMC`) but is a substantial computation: the
# original run needed a 900-row stratified subset because the original author found
# N > ~1500 rows already impractical for bayMDS's MCMC sampler on the full N=4256 (or,
# a fortiori, the production N=4586) player-seasons.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/archive/prototype_bmds_analysis.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c("dplyr", "ggplot2"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

RUN_BMDS_MCMC <- FALSE  # set TRUE only if you intend to re-run the (expensive) MCMC

dist_obj_path <- file.path(paths$processed, "archive_p23_distance_obj.rds")
if (!file.exists(dist_obj_path)) {
  stop("Run archive/exploratory_p23_similarity_analysis.R first to build ", dist_obj_path, call. = FALSE)
}
dist_obj <- readRDS(dist_obj_path)

meta <- dist_obj$row_keys |> mutate(season_id = as.factor(season_id))
D_use <- dist_obj$D_cosine

# ---------------------------------------------------------------------------
# 1) Stratified subset for BMDS (same construction as the original notebook)
# ---------------------------------------------------------------------------
min_minutes <- 900
N_target <- 900
meta0 <- meta |> mutate(row_id = row_number())
meta_f0 <- meta0 |> filter(minutes >= min_minutes)
n_seasons <- n_distinct(meta_f0$season_id)
n_per_season <- ceiling(N_target / n_seasons)

set.seed(69)
idx_keep <- meta_f0 |>
  group_by(season_id) |>
  group_modify(~ slice_sample(.x, n = min(nrow(.x), n_per_season))) |>
  ungroup()
keep_rows <- idx_keep$row_id

D_full <- as.matrix(D_use)
D_sub <- D_full[keep_rows, keep_rows]
meta_sub <- meta0[keep_rows, ] |> mutate(season_id = as.factor(season_id))

# ---------------------------------------------------------------------------
# 2) Classical MDS sanity check on the subset (cheap, always run)
# ---------------------------------------------------------------------------
mds2_sub <- cmdscale(D_sub, k = 2)
mds_sub_df <- meta_sub |> bind_cols(as.data.frame(mds2_sub) |> setNames(c("Dim1", "Dim2")))

p_cmds_sub <- ggplot(mds_sub_df, aes(Dim1, Dim2, color = season_id)) +
  geom_point(alpha = 0.6, size = 1.3) +
  theme_minimal() +
  labs(title = "Classical MDS (subset) on cosine dissimilarity")
ggsave(file.path(paths$figures, "archive_bmds_prototype_cmds_subset.png"), p_cmds_sub, width = 7, height = 5.5, dpi = 300)

# ---------------------------------------------------------------------------
# 3) Bayesian MDS: load the existing completed fit, or (optionally) re-run it
# ---------------------------------------------------------------------------
existing_fit_path <- file.path(dirname(dirname(paths$scripts)), "tfm_playground", "artifacts",
                                "bmds_cosine_subset_N900_p2_20260218_143220.rds")

if (RUN_BMDS_MCMC) {
  check_packages("bayMDS")
  library(bayMDS)
  set.seed(69)
  res_bmds <- bayMDS::bmdsMCMC(D_sub, p = 2, nwarm = 250, niter = 1000)
} else if (file.exists(existing_fit_path)) {
  b <- readRDS(existing_fit_path)
  res_bmds <- b$fit
  message("Loaded existing BMDS fit from ", existing_fit_path)
} else {
  stop("No existing BMDS fit found at ", existing_fit_path,
       " and RUN_BMDS_MCMC is FALSE. Set RUN_BMDS_MCMC <- TRUE to compute a fresh fit.",
       call. = FALSE)
}

coords_bmds <- as.data.frame(res_bmds$x_bmds)
colnames(coords_bmds) <- c("BMDS1", "BMDS2")
bmds_df <- meta_sub |> bind_cols(coords_bmds) |> mutate(r = sqrt(BMDS1^2 + BMDS2^2))

p_bmds <- ggplot(bmds_df, aes(BMDS1, BMDS2, color = season_id)) +
  geom_point(alpha = 0.65, size = 1.3) +
  theme_minimal() +
  labs(title = "BMDS p=2 (subset) on cosine dissimilarity")
ggsave(file.path(paths$figures, "archive_bmds_prototype_p2.png"), p_bmds, width = 7, height = 5.5, dpi = 300)

# Goodness-of-fit: observed vs. embedding-induced distances
D_hat <- as.matrix(dist(coords_bmds, method = "euclidean"))
get_upper <- function(M) M[upper.tri(M)]
cor_pearson <- cor(get_upper(D_sub), get_upper(D_hat), method = "pearson")
cor_spearman <- cor(get_upper(D_sub), get_upper(D_hat), method = "spearman")

writeLines(
  c(
    sprintf("N (subset) = %d", length(keep_rows)),
    sprintf("Radius (r) summary: mean=%.3f, sd=%.3f", mean(bmds_df$r), sd(bmds_df$r)),
    sprintf("Cor(d_obs, d_hat): Pearson=%.3f, Spearman=%.3f", cor_pearson, cor_spearman)
  ),
  file.path(paths$notes, "archive_bmds_prototype_summary.txt")
)
message("Prototype BMDS analysis complete. See outputs/notes/archive_bmds_prototype_summary.txt")
