# Stage 20c — fix stage 20b: stage 28's original CSVs have log_lik declared in the
# header but filled with placeholder zeros (compute_log_lik=0 at fit time, confirmed:
# stage 20b's K=3 loo() came back exactly zero/degenerate for every observation).
# Block C (scripts/C_block_c_loo_k2_vs_k3.R) already solved this by re-running
# generate_quantities() with compute_log_lik=1 against the stage 28 draws; that cached
# output (fits/csv/block_c_gq/k3/k3_gq-*.csv) has real, varying log_lik and is reused
# here directly rather than re-running generate_quantities.
#
# Diagonal (stage 03) log_lik was already real in its original CSVs (confirmed: stage
# 20b's diagonal loo() gave a plausible, varying, non-degenerate result) -- re-used
# via the same extraction as stage 20b.
#
# Overwrites: outputs/tables/20_icc_diagonal_vs_lowrank_rebuilt.csv
#             outputs/tables/20_pareto_k_diagonal_vs_lowrank_rebuilt.csv

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/20c_fix_diagonal_vs_k3_loo.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))
check_packages(c(required_base_packages, "loo"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(loo) })

py_script  <- file.path(model_root, "src", "extract_stan_csv_params.py")
SCRATCHPAD <- Sys.getenv("CLAUDE_SCRATCHPAD", unset = file.path(tempdir(), "20c_scratch"))
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

N <- 4586L

find_header <- function(csv_path) {
  con <- file(csv_path, "r"); on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0L) stop("EOF before header: ", csv_path)
    if (!startsWith(line, "#")) return(strsplit(line, ",", fixed = TRUE)[[1]])
  }
}

extract_log_lik_from <- function(csv_files, tag) {
  hdr <- find_header(csv_files[1])
  ll_names <- paste0("log_lik.", seq_len(N))
  idx0 <- match(ll_names, hdr) - 1L
  if (any(is.na(idx0))) stop("Missing log_lik columns for ", tag)
  names(idx0) <- ll_names

  chain_dfs <- lapply(seq_along(csv_files), function(i) {
    out <- file.path(SCRATCHPAD, sprintf("%s_chain%d.csv", tag, i))
    t0 <- proc.time()[["elapsed"]]
    args_pairs <- paste0(names(idx0), ":", idx0)
    ret <- system2("python3",
                   c(shQuote(py_script), shQuote(csv_files[i]), as.character(i),
                     shQuote(out), args_pairs),
                   stdout = FALSE, stderr = FALSE)
    if (ret != 0) stop("Python extraction failed for chain ", i, " of ", tag)
    message(sprintf("  [%s] chain %d done in %.0f s", tag, i, proc.time()[["elapsed"]] - t0))
    read.csv(out, check.names = FALSE)
  })
  draws_df <- do.call(rbind, chain_dfs)
  n_iter   <- max(draws_df[[".iteration"]])
  n_chains <- max(draws_df[[".chain"]])

  ll_arr <- array(NA_real_, dim = c(n_iter, n_chains, N))
  for (ch in seq_len(n_chains)) {
    sub <- draws_df[draws_df[[".chain"]] == ch, ll_names, drop = FALSE]
    ll_arr[, ch, ] <- as.matrix(sub)
  }
  ll_arr
}

message("Re-extracting log_lik for stage 03 (diagonal, rebuilt) -- unchanged from 20b...")
csv_diag <- discover_cmdstan_csv_files("03_real_diagonal_additive")
ll_diag  <- extract_log_lik_from(csv_diag, "diag")

message("\nExtracting log_lik for stage 28 (K=3) from the Block-C generate_quantities cache...")
csv_k3_gq <- sort(list.files(file.path(paths$fits, "csv", "block_c_gq", "k3"),
                             pattern = "[.]csv$", full.names = TRUE))
if (length(csv_k3_gq) == 0L) stop("No cached generate_quantities CSVs found for K=3 -- ",
                                   "run scripts/C_block_c_loo_k2_vs_k3.R first.")
message("  Found ", length(csv_k3_gq), " gq chain files: ", paste(basename(csv_k3_gq), collapse=", "))
ll_k3 <- extract_log_lik_from(csv_k3_gq, "k3gq")

message("\nComputing PSIS-LOO for both models...")
r_eff_diag <- relative_eff(exp(ll_diag), cores = 1)
loo_diag   <- loo(ll_diag, r_eff = r_eff_diag, cores = 1)

r_eff_k3 <- relative_eff(exp(ll_k3), cores = 1)
loo_k3   <- loo(ll_k3, r_eff = r_eff_k3, cores = 1)

message("\n=== Diagonal (stage 03, rebuilt) ===")
print(loo_diag)
message("\n=== K=3 Student-t (stage 28, rebuilt, via Block-C generate_quantities) ===")
print(loo_k3)

cmp <- loo_compare(list(diagonal_K0 = loo_diag, lowrank_K3_t = loo_k3))
message("\n=== loo_compare ===")
print(cmp)

cmp_df <- as.data.frame(cmp) |> tibble::rownames_to_column("model")
write_csv(cmp_df, file.path(paths$tables, "20_icc_diagonal_vs_lowrank_rebuilt.csv"))
message("\nSaved (corrected): outputs/tables/20_icc_diagonal_vs_lowrank_rebuilt.csv")

pareto_k_tbl <- tibble::tibble(
  model = c("diagonal_K0", "lowrank_K3_t"),
  frac_k_gt_0.5 = c(mean(loo_diag$diagnostics$pareto_k > 0.5),
                    mean(loo_k3$diagnostics$pareto_k > 0.5)),
  frac_k_gt_0.7 = c(mean(loo_diag$diagnostics$pareto_k > 0.7),
                    mean(loo_k3$diagnostics$pareto_k > 0.7))
)
write_csv(pareto_k_tbl, file.path(paths$tables, "20_pareto_k_diagonal_vs_lowrank_rebuilt.csv"))
message("Saved (corrected): outputs/tables/20_pareto_k_diagonal_vs_lowrank_rebuilt.csv")

message(sprintf("\nELPD diff (K3_t vs diagonal): %.1f (SE %.1f)",
                cmp["lowrank_K3_t", "elpd_diff"], cmp["lowrank_K3_t", "se_diff"]))
message("Stage 20c complete.")
