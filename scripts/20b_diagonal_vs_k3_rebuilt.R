# Stage 20b — LOO comparison, diagonal (K=0) vs low-rank Student-t (K=3), on the
# REBUILT dataset (P=48, N=4586).
#
# The original scripts/20_diagonal_comparison.R compares stage 03 vs stage 10
# (K=2, Normal) and was never rebuilt after the P=48 rebuild; per CLAUDE.md
# Issue 7 / "Next required go-ahead", the wanted comparison is the actual
# PRODUCTION model (stage 28, K=3, Student-t) against the diagonal baseline,
# both refit on rebuilt data (both already have log_lik computed -- confirmed
# by inspecting their CSV headers directly, so no new Stan fit is needed here).
#
# log_lik is extracted via the same fast Python column-pass approach used
# elsewhere in this project (src/extract_stan_csv_params.py), rather than
# cmdstanr::as_cmdstan_fit(), to avoid loading the full ~79K/~236K-column CSVs.
#
# Writes: outputs/tables/20_icc_diagonal_vs_lowrank_rebuilt.csv (loo_compare summary)
#         outputs/tables/20_pareto_k_diagonal_vs_lowrank_rebuilt.csv (Pareto-k diagnostics)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/20b_diagonal_vs_k3_rebuilt.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))
check_packages(c(required_base_packages, "loo"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(loo) })

py_script  <- file.path(model_root, "src", "extract_stan_csv_params.py")
SCRATCHPAD <- Sys.getenv("CLAUDE_SCRATCHPAD", unset = file.path(tempdir(), "20b_scratch"))
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

extract_log_lik <- function(model_id) {
  csv_files <- discover_cmdstan_csv_files(model_id)
  if (length(csv_files) == 0L) stop("No CSV chains for ", model_id)
  message(sprintf("[%s] %d chain(s): %s", model_id, length(csv_files),
                   paste(basename(dirname(csv_files[1])))))

  hdr <- find_header(csv_files[1])
  ll_names <- paste0("log_lik.", seq_len(N))
  idx0 <- match(ll_names, hdr) - 1L
  if (any(is.na(idx0))) stop("Missing log_lik columns for ", model_id)
  names(idx0) <- ll_names

  chain_dfs <- lapply(seq_along(csv_files), function(i) {
    out <- file.path(SCRATCHPAD, sprintf("%s_chain%d.csv", model_id, i))
    t0 <- proc.time()[["elapsed"]]
    args_pairs <- paste0(names(idx0), ":", idx0)
    ret <- system2("python3",
                   c(shQuote(py_script), shQuote(csv_files[i]), as.character(i),
                     shQuote(out), args_pairs),
                   stdout = FALSE, stderr = FALSE)
    if (ret != 0) stop("Python extraction failed for chain ", i, " of ", model_id)
    message(sprintf("  chain %d done in %.0f s", i, proc.time()[["elapsed"]] - t0))
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

message("Extracting log_lik for stage 03 (diagonal, rebuilt)...")
t0 <- proc.time()[["elapsed"]]
ll_diag <- extract_log_lik("03_real_diagonal_additive")
message(sprintf("Total: %.0f s", proc.time()[["elapsed"]] - t0))

message("\nExtracting log_lik for stage 28 (K=3, Student-t, rebuilt)...")
t0 <- proc.time()[["elapsed"]]
ll_k3 <- extract_log_lik("28_real_lowrank_a_diag_b_t_k3")
message(sprintf("Total: %.0f s", proc.time()[["elapsed"]] - t0))

message("\nComputing PSIS-LOO for both models...")
r_eff_diag <- relative_eff(exp(ll_diag), cores = 1)
loo_diag   <- loo(ll_diag, r_eff = r_eff_diag, cores = 1)

r_eff_k3 <- relative_eff(exp(ll_k3), cores = 1)
loo_k3   <- loo(ll_k3, r_eff = r_eff_k3, cores = 1)

message("\n=== Diagonal (stage 03, rebuilt) ===")
print(loo_diag)
message("\n=== K=3 Student-t (stage 28, rebuilt) ===")
print(loo_k3)

cmp <- loo_compare(list(diagonal_K0 = loo_diag, lowrank_K3_t = loo_k3))
message("\n=== loo_compare ===")
print(cmp)

cmp_df <- as.data.frame(cmp) |> tibble::rownames_to_column("model")
write_csv(cmp_df, file.path(paths$tables, "20_icc_diagonal_vs_lowrank_rebuilt.csv"))
message("\nSaved: outputs/tables/20_icc_diagonal_vs_lowrank_rebuilt.csv")

pareto_k_tbl <- tibble::tibble(
  model = c("diagonal_K0", "lowrank_K3_t"),
  frac_k_gt_0.5 = c(mean(loo_diag$diagnostics$pareto_k > 0.5),
                    mean(loo_k3$diagnostics$pareto_k > 0.5)),
  frac_k_gt_0.7 = c(mean(loo_diag$diagnostics$pareto_k > 0.7),
                    mean(loo_k3$diagnostics$pareto_k > 0.7))
)
write_csv(pareto_k_tbl, file.path(paths$tables, "20_pareto_k_diagonal_vs_lowrank_rebuilt.csv"))
message("Saved: outputs/tables/20_pareto_k_diagonal_vs_lowrank_rebuilt.csv")

message(sprintf("\nELPD diff (K3_t vs diagonal): %.1f (SE %.1f)",
                cmp["lowrank_K3_t", "elpd_diff"], cmp["lowrank_K3_t", "se_diff"]))
message("Stage 20b complete.")
