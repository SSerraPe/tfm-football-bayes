# B4 post-fit diagnostics: Task L (did nu move?) + Task S (did Lambda_a mixing improve?)
#
# Run immediately after B3 (stages 10 + 18 rebuild) completes.
# Compares against pre-rebuild baselines:
#   nu: 2.95
#   Lambda_a[10,1] raw ESS=19, Rhat=1.185
#   LLt.1.10 ESS=24, Rhat=1.130
#   Post-GK-only ablation baseline: LLt.1.10 ESS=64, Rhat=1.048

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/B4_postfit_diagnostics.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "posterior"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(posterior) })

py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")

SCRATCHPAD <- file.path(
  "/private/tmp/claude-501/-Users-sserra-Documents-01-MESIO-UPC-TFM-model",
  "76818d40-38b9-4367-9dc1-5d60a0b70149", "scratchpad"
)
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

# ---- Find stage-18 rebuilt CSV files ------------------------------------------

MODEL_ID <- "18_real_lowrank_a_diag_b_t"
csv_files <- discover_cmdstan_csv_files(MODEL_ID)
if (length(csv_files) == 0) stop("No CSV files found for ", MODEL_ID, call. = FALSE)
message("Found ", length(csv_files), " chain CSVs for ", MODEL_ID)
message("  Latest: ", csv_files[1])

# ---- Column finder (same as in stage 26) -------------------------------------

parse_col_index <- function(csv_path, param_name) {
  con <- file(csv_path, "r")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0L) return(NA_integer_)
    if (!startsWith(line, "#")) {
      hdr <- strsplit(line, ",")[[1]]
      idx <- which(hdr == param_name)
      return(if (length(idx) == 1L) idx[1L] else NA_integer_)
    }
  }
}

extract_param_chains <- function(csv_files, param_name, scratchpad = SCRATCHPAD) {
  col_r <- parse_col_index(csv_files[1], param_name)
  if (is.na(col_r)) stop("Parameter '", param_name, "' not found in CSV header", call. = FALSE)
  message("  ", param_name, " is at column ", col_r, " (1-indexed)")
  col_0 <- col_r - 1L
  chain_dfs <- lapply(seq_along(csv_files), function(i) {
    tmp <- file.path(scratchpad, sprintf("b4_%d_%s.csv", i, gsub("\\.", "_", param_name)))
    system2("python3", args = c(shQuote(py_script), shQuote(csv_files[i]),
                                as.character(i), shQuote(tmp),
                                paste0(param_name, ":", col_0)),
            stdout = FALSE, stderr = FALSE)
    read.csv(tmp, check.names = FALSE)
  })
  do.call(rbind, chain_dfs)
}

# ---- Task L: nu ---------------------------------------------------------------

message("\n=== Task L: Did nu move? (pre-rebuild baseline: nu ~ 2.95) ===")

df_nu <- extract_param_chains(csv_files, "nu")
chains_nu <- df_nu[[".chain"]]
nu_vals   <- df_nu[["nu"]]

arr_nu <- {
  n_iter <- max(df_nu[[".iteration"]])
  a <- array(NA_real_, dim = c(n_iter, 4L, 1L))
  for (ch in 1:4) a[, ch, 1] <- nu_vals[chains_nu == ch]
  as_draws_array(a)
}

nu_mean <- mean(nu_vals)
nu_ess  <- ess_bulk(arr_nu)
nu_rhat <- rhat(arr_nu)

message(sprintf("  nu grand mean = %.3f  (pre-rebuild: 2.95)", nu_mean))
message(sprintf("  nu ESS_bulk   = %.0f   Rhat = %.3f", nu_ess, nu_rhat))

task_l_verdict <- if (nu_mean >= 4.5) {
  "TRANSFORMS ABSORBED HEAVY TAILS (nu >= 4.5) — t-model still appropriate."
} else {
  "NU STILL PINNED NEAR 3 (nu < 4.5) — two-component mixture noted for future work; proceeding."
}
message("  Task L verdict: ", task_l_verdict)

# ---- Task S: Lambda_a[10,1] and LLt.1.10 mixing ------------------------------

message("\n=== Task S: Did Lambda_a mixing improve? ===")
message("  Pre-rebuild baselines: Lambda_a ESS=19/Rhat=1.185; LLt.1.10 ESS=24/Rhat=1.130")
message("  Post-GK-only ablation: LLt.1.10 ESS=64/Rhat=1.048")

df_la10 <- extract_param_chains(csv_files, "Lambda_a.10.1")
df_la1  <- extract_param_chains(csv_files, "Lambda_a.1.1")

chains_la <- df_la10[[".chain"]]
la10_1    <- df_la10[["Lambda_a.10.1"]]
la1_1     <- df_la1[["Lambda_a.1.1"]]
llt_1_10  <- la1_1 * la10_1

mk_arr <- function(x, chains) {
  n_iter <- max(df_la10[[".iteration"]])
  a <- array(NA_real_, dim = c(n_iter, 4L, 1L))
  for (ch in 1:4) a[, ch, 1] <- x[chains == ch]
  as_draws_array(a)
}

arr_la10 <- mk_arr(la10_1,   chains_la)
arr_llt  <- mk_arr(llt_1_10, chains_la)

ess_la10  <- ess_bulk(arr_la10);  rht_la10  <- rhat(arr_la10)
ess_llt   <- ess_bulk(arr_llt);   rht_llt   <- rhat(arr_llt)

chain_means_10_1 <- tapply(la10_1, chains_la, mean)
between_sd <- sd(chain_means_10_1)

message(sprintf("  Lambda_a[10,1] grand mean = %.4f  between-chain SD = %.4f",
                mean(la10_1), between_sd))
message(sprintf("  Lambda_a[10,1] ESS=%.0f  Rhat=%.3f  (pre-rebuild ESS=19, Rhat=1.185)",
                ess_la10, rht_la10))
message(sprintf("  LLt.1.10       ESS=%.0f  Rhat=%.3f  (pre-rebuild ESS=24, Rhat=1.130)",
                ess_llt, rht_llt))

task_s_verdict <- if (ess_llt > 200 && rht_llt < 1.05) {
  "MIXING RESOLVED (ESS > 200 and Rhat < 1.05) — proceed to Task O."
} else {
  paste0("MIXING STILL POOR (ESS=", round(ess_llt), ", Rhat=", round(rht_llt, 3),
         ") — escalated; proceed to Task O with flag.")
}
message("  Task S verdict: ", task_s_verdict)

# ---- Save summary ------------------------------------------------------------

result <- tibble(
  task          = c("L_nu", "S_LLt"),
  metric        = c(sprintf("nu=%.3f", nu_mean),
                    sprintf("ESS=%.0f,Rhat=%.3f", ess_llt, rht_llt)),
  baseline      = c("nu=2.95", "ESS=24,Rhat=1.130"),
  verdict       = c(task_l_verdict, task_s_verdict)
)
write_csv(result, file.path(paths$tables, "B4_task_L_S_verdicts.csv"))
message("\nSaved B4_task_L_S_verdicts.csv")
message("\nB4 diagnostics complete.")
