# Stage 26 — GK-only ablation fit: Student-t model with GK rows excluded.
#
# Purpose: isolate the effect of GK exclusion on Lambda_a mixing, specifically
# Lambda_a[10,1] (per90_passes loading on factor 1), which drove the poor ESS
# in the LLt.1.10 off-diagonal entry (pre-GK ESS=24, Rhat=1.13).
#
# This script does NOT touch any production files (model_objects.rds, Y_scaled.csv,
# stan_data/). It builds a GK-excluded dataset on the fly and fits with a distinct
# model_id "26_ablation_gk_only_t".
#
# Pre-GK baseline (from stage-18 draws, per chain):
#   Chain 1: Lambda_a[10,1] mean = -0.4476
#   Chain 2: Lambda_a[10,1] mean = -0.4325
#   Chain 3: Lambda_a[10,1] mean = -0.4627
#   Chain 4: Lambda_a[10,1] mean = -0.4657
#   Between-chain SD = 0.0153  (ESS of LLt.1.10 = 24, Rhat = 1.13)
#
# Prediction: if GK bimodality drives the mixing failure, removing 358 GK rows
# / 121 GK players will concentrate Lambda_a[10,1]'s posterior and the
# four chains' means will converge tightly.
#
# Post-fit diagnostic (run immediately after sampling):
#   - Per-chain Lambda_a[10,1] means (primary endpoint)
#   - ESS / Rhat for Lambda_a[10,1] and LLt.1.10 (rotation-invariant cross-product)
#
# Output:
#   fits/csv/26_ablation_gk_only_t/<timestamp>/*.csv
#   outputs/tables/26_ablation_lambda_chain_means.csv
#   outputs/notes/note_ablation_gk_lambda.md

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/26_ablation_gk_only_t.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr", "posterior"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(posterior)
})

MODEL_ID  <- "26_real_lowrank_a_diag_b_t_gk_excl"
RANK_A    <- 2L  # same as stage 18

# ---- 1. Load production model objects and Y_scaled ---------------------------

mo        <- readRDS(file.path(paths$processed, "model_objects.rds"))
Y_full    <- as.matrix(read_csv(file.path(paths$processed, "Y_scaled.csv"),
                                 show_col_types = FALSE))
long_raw  <- read_csv(file.path(paths$processed, "longitudinal_player_seasons.csv"),
                      col_types = cols_only(
                        player_id = col_double(),
                        season_id = col_double(),
                        total_gk_saves = col_double()
                      ))

stopifnot(nrow(Y_full) == nrow(mo$metadata))
message("Production dataset: N=", nrow(Y_full), ", I=", length(mo$player_levels),
        ", P=", ncol(Y_full))

# ---- 2. Build GK flag at player-season level (total_gk_saves > 0) -----------

gk_flag <- long_raw |>
  mutate(is_gk = !is.na(total_gk_saves) & total_gk_saves > 0) |>
  select(player_id, season_id, is_gk)

meta_flagged <- mo$metadata |>
  left_join(gk_flag, by = c("player_id", "season_id")) |>
  mutate(is_gk = coalesce(is_gk, FALSE))

n_gk_rows  <- sum(meta_flagged$is_gk)
n_gk_plrs  <- n_distinct(meta_flagged$player_id[meta_flagged$is_gk])
message("GK rows to remove: ", n_gk_rows, " (", n_gk_plrs, " distinct GK players)")

# ---- 3. Filter rows and remap indices ----------------------------------------

keep      <- !meta_flagged$is_gk
Y_abl     <- Y_full[keep, ]
meta_abl  <- meta_flagged[keep, ]

# Remap player_index: keep only players who still have at least one row
# mo$player_levels[j] gives the player_id corresponding to old index j
old_pid_by_row   <- mo$player_levels[mo$player_index]   # player_id for each row
pid_keep_rows    <- old_pid_by_row[keep]                 # player_ids for kept rows
new_player_levels <- unique(pid_keep_rows)               # in occurrence order
# Use mo$player_levels ordering for stability
new_player_levels <- mo$player_levels[mo$player_levels %in% new_player_levels]

new_player_index  <- match(pid_keep_rows, new_player_levels)
new_season_index  <- mo$season_index[keep]               # seasons unchanged

stopifnot(all(!is.na(new_player_index)))
stopifnot(nrow(Y_abl) == sum(keep))

message("Ablation dataset: N=", nrow(Y_abl), ", I=", length(new_player_levels),
        ", S=", length(mo$season_levels), ", P=", ncol(Y_abl))

# ---- 4. Build ablation model_objects and stan_data --------------------------

mo_abl <- mo
mo_abl$metadata      <- meta_abl
mo_abl$player_index  <- new_player_index
mo_abl$season_index  <- new_season_index
mo_abl$player_levels <- new_player_levels

stan_data_abl <- build_lowrank_a_diag_b_stan_data(
  model_objects = mo_abl,
  Y             = Y_abl,
  rank_a        = RANK_A
)
# Ablation: skip log_lik to keep CSVs smaller (we don't need LOO for this fit)
stan_data_abl$compute_log_lik <- 0L
message("Stan data built. N=", stan_data_abl$N, " I=", stan_data_abl$I,
        " P=", stan_data_abl$P)

# ---- 5. PCA init (same as stage 18) -----------------------------------------

base_init_fn <- build_pca_init(stan_data_abl)
init_fn <- function(chain_id = 1) {
  init <- base_init_fn(chain_id)
  init$mu <- NULL
  init$nu <- 30.0
  init
}

# ---- 6. Fit ---------------------------------------------------------------

fit_path <- file.path(paths$fits, paste0(MODEL_ID, "_fit.rds"))

fit <- fit_cmdstan_model(
  stan_file  = file.path(paths$stan, "additive_lowrank_a_diag_b_t.stan"),
  stan_data  = stan_data_abl,
  fit_path   = fit_path,
  seed       = pipeline$seed + 260L,
  model_id   = MODEL_ID,
  init       = init_fn
)

if (is.null(fit)) {
  message("Sampling skipped (BFA_RUN_STAN=false). Set BFA_RUN_STAN=true to run.")
  quit(save = "no", status = 0)
}

write_fit_outputs(fit, MODEL_ID)

# ---- 7. Post-fit diagnostic: Lambda_a[10,1] per-chain means -----------------
# Parse CSV header to find column index for Lambda_a.10.1
# CmdStan CSVs: first ~47 comment lines (#), then the header row.

csv_files <- fit$output_files()

parse_col_index <- function(csv_path, param_name) {
  con <- file(csv_path, "r")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) return(NA_integer_)
    if (!startsWith(line, "#")) {
      header <- strsplit(line, ",")[[1]]
      idx <- which(header == param_name)
      return(if (length(idx) == 1L) idx[1L] else NA_integer_)
    }
  }
}

# R (1-indexed) position of Lambda_a.10.1
col_r <- parse_col_index(csv_files[1], "Lambda_a.10.1")
col_r_la11 <- parse_col_index(csv_files[1], "Lambda_a.1.1")
message("Lambda_a.10.1 is at column ", col_r, " (1-indexed) in the CSV")
message("Lambda_a.1.1  is at column ", col_r_la11, " (1-indexed) in the CSV")

# Python column-pass for Lambda_a.10.1 and Lambda_a.1.1 across all chains
py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")

extract_param_chains <- function(csv_files, col_r, param_name, scratchpad) {
  col_0idx <- col_r - 1L   # Python uses 0-indexed
  chain_dfs <- lapply(seq_along(csv_files), function(i) {
    tmp <- file.path(scratchpad, sprintf("abl_chain%d_%s.csv", i, gsub("\\.", "_", param_name)))
    system2("python3",
      args = c(shQuote(py_script), shQuote(csv_files[i]), as.character(i),
               shQuote(tmp), paste0(param_name, ":", col_0idx)),
      stdout = FALSE, stderr = FALSE)
    df <- read.csv(tmp, check.names = FALSE)
    df
  })
  do.call(rbind, chain_dfs)
}

scratchpad_dir <- Sys.getenv("CLAUDE_SCRATCHPAD",
  "/private/tmp/claude-501/-Users-sserra-Documents-01-MESIO-UPC-TFM-model/76818d40-38b9-4367-9dc1-5d60a0b70149/scratchpad")
dir.create(scratchpad_dir, recursive = TRUE, showWarnings = FALSE)

message("\nExtracting Lambda_a.10.1 draws across chains...")
df_la10_1 <- extract_param_chains(csv_files, col_r, "Lambda_a.10.1", scratchpad_dir)
df_la1_1  <- extract_param_chains(csv_files, col_r_la11, "Lambda_a.1.1", scratchpad_dir)

chains <- df_la10_1[[".chain"]]
la10_1 <- df_la10_1[["Lambda_a.10.1"]]
la1_1  <- df_la1_1[["Lambda_a.1.1"]]

# Per-chain means
chain_means_10_1 <- tapply(la10_1, chains, mean)
chain_sds_10_1   <- tapply(la10_1, chains, sd)
chain_means_1_1  <- tapply(la1_1, chains, mean)

# LLt.1.10 per draw
llt_1_10 <- la1_1 * la10_1

# ESS / Rhat for Lambda_a[10,1] and LLt.1.10
mk_arr <- function(x, chains, n_iter = 1000L) {
  arr <- array(NA_real_, dim = c(n_iter, 4L, 1L))
  for (ch in 1:4) arr[, ch, 1] <- x[chains == ch]
  arr
}
n_iter <- max(df_la10_1[[".iteration"]])
arr_la10 <- as_draws_array(mk_arr(la10_1, chains, n_iter))
arr_llt  <- as_draws_array(mk_arr(llt_1_10, chains, n_iter))

ess_la10 <- ess_bulk(arr_la10)
rht_la10 <- rhat(arr_la10)
ess_llt  <- ess_bulk(arr_llt)
rht_llt  <- rhat(arr_llt)

message("\n=== Post-GK ablation: Lambda_a[10,1] (per90_passes, factor 1) ===")
for (ch in 1:4) {
  message(sprintf("  Chain %d: mean=%+.4f  sd=%.4f", ch, chain_means_10_1[ch], chain_sds_10_1[ch]))
}
message(sprintf("  Grand mean: %+.4f | between-chain SD: %.4f",
        mean(la10_1), sd(chain_means_10_1)))
message(sprintf("  ESS_bulk = %.0f | Rhat = %.3f", ess_la10, rht_la10))

message(sprintf("\nLLt.1.10 ESS_bulk = %.0f | Rhat = %.3f  (pre-GK: ESS=24, Rhat=1.130)",
        ess_llt, rht_llt))

# ---- 8. Save results --------------------------------------------------------

pre_gk_means <- c(-0.4476, -0.4325, -0.4627, -0.4657)  # from scratchpad pre-GK analysis

result_tbl <- data.frame(
  chain        = 1:4,
  pre_gk_mean  = pre_gk_means,
  post_gk_mean = as.numeric(chain_means_10_1),
  change       = as.numeric(chain_means_10_1) - pre_gk_means,
  post_gk_sd   = as.numeric(chain_sds_10_1)
)

write.csv(result_tbl,
  file.path(paths$tables, "26_ablation_lambda_chain_means.csv"), row.names = FALSE)

# Verdict
between_sd_pre  <- 0.0153
between_sd_post <- sd(chain_means_10_1)
convergence_ratio <- between_sd_post / between_sd_pre

note_lines <- c(
  "# Ablation Note: GK-Only Exclusion Effect on Lambda_a[10,1]",
  "",
  "## Pre-GK baseline (stage-18 production fit)",
  "",
  "| Chain | Lambda_a[10,1] mean | between-chain SD |",
  "|---|---|---|",
  sprintf("| 1 | %.4f | — |", pre_gk_means[1]),
  sprintf("| 2 | %.4f | — |", pre_gk_means[2]),
  sprintf("| 3 | %.4f | — |", pre_gk_means[3]),
  sprintf("| 4 | %.4f | — |", pre_gk_means[4]),
  sprintf("| **combined** | **%.4f** | **%.4f** |", mean(pre_gk_means), between_sd_pre),
  "",
  "LLt.1.10 ESS=24, Rhat=1.130.",
  "",
  "## Post-GK ablation (stage-26, GK rows excluded)",
  "",
  sprintf("N reduced: 4944 → %d (removed %d rows / %d GK players)",
          nrow(Y_abl), n_gk_rows, n_gk_plrs),
  sprintf("I reduced: 1650 → %d", length(new_player_levels)),
  "",
  "| Chain | Lambda_a[10,1] mean | between-chain SD |",
  "|---|---|---|",
  sprintf("| 1 | %.4f | — |", chain_means_10_1[1]),
  sprintf("| 2 | %.4f | — |", chain_means_10_1[2]),
  sprintf("| 3 | %.4f | — |", chain_means_10_1[3]),
  sprintf("| 4 | %.4f | — |", chain_means_10_1[4]),
  sprintf("| **combined** | **%.4f** | **%.4f** |", mean(la10_1), between_sd_post),
  "",
  sprintf("LLt.1.10 ESS=%.0f, Rhat=%.3f.", ess_llt, rht_llt),
  sprintf("Lambda_a[10,1] ESS=%.0f, Rhat=%.3f.", ess_la10, rht_la10),
  "",
  "## Verdict",
  "",
  sprintf("Between-chain SD: %.4f → %.4f (ratio %.2f×).",
          between_sd_pre, between_sd_post, convergence_ratio),
  if (convergence_ratio < 0.4) {
    "**GK-driven bimodality confirmed**: chains converge substantially after GK exclusion. Proceed with Task K full bundle."
  } else if (convergence_ratio < 0.7) {
    "**Partial GK effect**: some improvement but mixing issue not fully resolved by GK exclusion alone. Full bundle (transforms) may help further; investigate which features drive residual disagreement."
  } else {
    "**GK exclusion did not resolve mixing**: between-chain disagreement unchanged. The per90_passes mixing problem has a different cause. Do not attribute to GK bimodality."
  },
  "",
  sprintf("_Written by scripts/26_ablation_gk_only_t.R on %s_", Sys.Date())
)

writeLines(note_lines, file.path(paths$notes, "note_ablation_gk_lambda.md"))
message("\nSaved: 26_ablation_lambda_chain_means.csv, note_ablation_gk_lambda.md")
message("\nStage 26 complete.")
