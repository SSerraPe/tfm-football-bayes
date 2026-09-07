# Stage 21b — Build the missing full posterior summary for stage 21 (K=4).
#
# stage 21's fit completed (4 chains of real posterior draws exist under
# fits/csv/21_real_k4_t/...) but outputs/tables/21_real_k4_t_posterior_summary.csv
# only ever got a 1-row (lp__-only) summary written -- the full per-parameter
# extraction (nu, Lambda_a, psi_a, sigma_b, sigma_e) that scripts/29b and 31b
# require never ran. This reproduces it via the same Python column-pass approach
# used for stage 28 (scripts/B_block_b_diagnostics_28.R), generalised to K=4.
#
# Writes: outputs/tables/21_real_k4_t_posterior_summary.csv
#         (same variable/mean/median/sd/mad/q5/q95/rhat/ess_bulk/ess_tail format
#          as the existing stage 18/28 posterior summaries)

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else
  "model/scripts/21b_build_k4_posterior_summary.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))),
                 "src", "bootstrap.R"))
check_packages(c(required_base_packages, "posterior"))

suppressPackageStartupMessages({ library(readr); library(dplyr); library(posterior) })

py_script <- file.path(model_root, "src", "extract_stan_csv_params.py")
SCRATCHPAD <- Sys.getenv("CLAUDE_SCRATCHPAD",
                         unset = file.path(tempdir(), "21b_scratch"))
dir.create(SCRATCHPAD, recursive = TRUE, showWarnings = FALSE)

MODEL_ID <- "21_real_k4_t"
K <- 4L; P <- 48L

csv_files <- discover_cmdstan_csv_files(MODEL_ID)
if (length(csv_files) == 0) stop("No CSV files found for ", MODEL_ID)
message("Found ", length(csv_files), " chain CSVs for ", MODEL_ID)

find_col_indices <- function(csv_path, param_names) {
  con <- file(csv_path, "r")
  on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0L) return(setNames(rep(NA_integer_, length(param_names)), param_names))
    if (!startsWith(line, "#")) {
      hdr <- strsplit(line, ",")[[1]]
      idx <- match(param_names, hdr)
      return(setNames(idx, param_names))
    }
  }
}

message("Parsing CSV header to find column positions...")
lambda_names  <- paste0("Lambda_a.", rep(1:P, K), ".", rep(1:K, each = P))
psi_names     <- paste0("psi_a.", 1:P)
sigma_b_names <- paste0("sigma_b.", 1:P)
sigma_e_names <- paste0("sigma_e.", 1:P)
all_params    <- c("nu", lambda_names, psi_names, sigma_b_names, sigma_e_names)

col_idx <- find_col_indices(csv_files[1], all_params)
missing <- names(col_idx)[is.na(col_idx)]
if (length(missing) > 0) message("WARNING: not found in header: ", paste(missing, collapse=", "))
found_params <- col_idx[!is.na(col_idx)]
message("Found ", length(found_params), " of ", length(all_params), " parameters")

extract_batch <- function(csv_path, chain_id, params_0idx, out_path) {
  args_pairs <- paste0(names(params_0idx), ":", params_0idx - 1L)
  cmd_args <- c(shQuote(py_script), shQuote(csv_path),
                as.character(chain_id), shQuote(out_path),
                args_pairs)
  ret <- system2("python3", args = cmd_args, stdout = FALSE, stderr = FALSE)
  if (ret != 0) stop("Python extraction failed for chain ", chain_id)
  read.csv(out_path, check.names = FALSE)
}

message("Extracting ", length(found_params), " params from ", length(csv_files), " chains...")
t_all0 <- proc.time()[["elapsed"]]
chain_dfs <- lapply(seq_along(csv_files), function(i) {
  out <- file.path(SCRATCHPAD, sprintf("k4_chain%d.csv", i))
  t0 <- proc.time()[["elapsed"]]
  message("  Chain ", i, "...")
  df <- extract_batch(csv_files[i], i, found_params, out)
  message(sprintf("    done in %.0f s", proc.time()[["elapsed"]] - t0))
  df
})
message(sprintf("Total extraction time: %.0f s", proc.time()[["elapsed"]] - t_all0))

draws_df <- do.call(rbind, chain_dfs)
n_iter   <- max(draws_df[[".iteration"]])
n_chains <- max(draws_df[[".chain"]])
message(sprintf("Combined: %d iterations x %d chains = %d total draws",
                n_iter, n_chains, nrow(draws_df)))

param_cols <- setdiff(names(draws_df), c(".chain", ".iteration"))

make_array <- function(df, param, n_iter, n_chains) {
  a <- array(NA_real_, dim = c(n_iter, n_chains, 1L))
  for (ch in seq_len(n_chains)) a[, ch, 1L] <- df[[param]][df[[".chain"]] == ch]
  as_draws_array(a)
}

message("Computing summary stats for ", length(param_cols), " parameters...")
summ_rows <- lapply(param_cols, function(p) {
  v   <- draws_df[[p]]
  arr <- make_array(draws_df, p, n_iter, n_chains)
  data.frame(
    variable = p,
    mean     = mean(v),
    median   = median(v),
    sd       = sd(v),
    mad      = mad(v),
    q5       = quantile(v, 0.05, names = FALSE),
    q95      = quantile(v, 0.95, names = FALSE),
    rhat     = rhat(arr),
    ess_bulk = ess_bulk(arr),
    ess_tail = ess_tail(arr),
    row.names = NULL
  )
})
summ_df <- do.call(rbind, summ_rows)

# Rename dot-notation variables to the bracket notation used by
# posterior::summarise_draws() / the existing stage 18/28 summary files, so
# scripts/29b_kcompare_residuals.R and scripts/31b_kcompare_pca.R can parse them.
rename_var <- function(x) {
  if (x == "nu") return("nu")
  m2 <- regmatches(x, regexec("^([a-zA-Z_]+)\\.(\\d+)\\.(\\d+)$", x))[[1]]
  if (length(m2) == 4) return(sprintf("%s[%s,%s]", m2[2], m2[3], m2[4]))
  m1 <- regmatches(x, regexec("^([a-zA-Z_]+)\\.(\\d+)$", x))[[1]]
  if (length(m1) == 3) return(sprintf("%s[%s]", m1[2], m1[3]))
  x
}
summ_df$variable <- vapply(summ_df$variable, rename_var, character(1))

out_path <- file.path(paths$tables, "21_real_k4_t_posterior_summary.csv")
write_csv(summ_df, out_path)
message("\nSaved full posterior summary: ", out_path, " (", nrow(summ_df), " rows)")

nu_row <- summ_df[summ_df$variable == "nu", ]
message(sprintf("nu: mean=%.4f  95%%CI=[%.4f, %.4f]  ESS_bulk=%.0f  Rhat=%.4f",
                nu_row$mean, nu_row$q5, nu_row$q95, nu_row$ess_bulk, nu_row$rhat))
message("Stage 21b complete.")
