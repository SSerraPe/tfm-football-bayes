# Utilities for extracting parameters from large CmdStan CSV files without
# loading them via cmdstanr::as_cmdstan_fit(), which parses all columns and
# requires ~45 min and ~10 GB RAM for the stage-18 t-model (12 GB, 276,606 cols).
#
# Primary function: extract_stan_csv_params()
# Wraps src/extract_stan_csv_params.py via system2("python3", ...).
#
# Column index convention:
#   All column positions are 0-INDEXED, matching Python / the CSV header.
#   "col 8" in the project memory document means 1-indexed → use 7 here.
#   lp__ is index 0; sampler-internal cols are indices 1-6; first model
#   parameter starts at index 7.
#
# Documented positions for stage-18 PRE-REBUILD (additive_lowrank_a_diag_b_t, K=2, P=53):
# STALE after Session 8 rebuild (stage-18 now uses P=48). STAGE18_* constants below
# are preserved for reference but NOT USED anywhere. Use parse_col_index() or the
# dynamic discovery in B_block_b_diagnostics_28.R for any new extraction.
#   nu               : col 7
#   lambda_a_diag[k] : cols 90758-90759  (K=2)
#   psi_a[p]         : cols 90863-90915  (P=53)
#   sigma_b[p]       : cols 91552-91604
#   sigma_e[p]       : cols 91605-91657
#   Lambda_a[p,q]    : cols 91658-91763  (P*K=106, column-major)
#   Sigma_a_diag[p]  : cols 276288-276340
#   var_b[p]         : cols 276341-276393
#   var_e[p]         : cols 276394-276446
#   prop_a[p]        : cols 276447-276499
#   prop_b[p]        : cols 276500-276552
#   prop_e[p]        : cols 276553-276605
#
# To re-derive positions after a new fit (different P or K):
#   head_line <- readLines(csv_path, n = 55)
#   header    <- head_line[!startsWith(head_line, "#")][1]
#   cols      <- strsplit(header, ",")[[1]]
#   which(cols == "nu") - 1L  # subtract 1 for 0-indexed

# ---- Helper builders -------------------------------------------------------

#' Build a named integer vector for a contiguous range of 1-D array elements.
build_col_range <- function(name_prefix, start_col, n) {
  nms  <- paste0(name_prefix, ".", seq_len(n))
  cols <- as.integer(start_col) + seq.int(0L, n - 1L)
  structure(cols, names = nms)
}

#' Build a named integer vector for a 2-D P×K matrix (column-major, CmdStan).
build_col_range_2d <- function(name_prefix, start_col, P, K) {
  nms  <- paste0(name_prefix, ".",
                 rep(seq_len(P), times = K), ".",
                 rep(seq_len(K), each = P))
  cols <- as.integer(start_col) + seq.int(0L, P * K - 1L)
  structure(cols, names = nms)
}

# ---- Pre-built column maps for stage-18 ------------------------------------

STAGE18_NU_MAP      <- c(nu = 7L)
STAGE18_SIGMA_E_MAP <- build_col_range("sigma_e",  91605L, 53L)
STAGE18_SIGMA_B_MAP <- build_col_range("sigma_b",  91552L, 53L)
STAGE18_PSI_A_MAP   <- build_col_range("psi_a",    90863L, 53L)
STAGE18_LAMBDA_A_MAP <- build_col_range_2d("Lambda_a", 91658L, 53L, 2L)
STAGE18_LAMBDA_A_DIAG_MAP <- c(lambda_a_diag.1 = 90758L, lambda_a_diag.2 = 90759L)
STAGE18_PROP_A_MAP  <- build_col_range("prop_a",   276447L, 53L)
STAGE18_PROP_B_MAP  <- build_col_range("prop_b",   276500L, 53L)
STAGE18_PROP_E_MAP  <- build_col_range("prop_e",   276553L, 53L)
STAGE18_SIGMA_A_DIAG_MAP <- build_col_range("Sigma_a_diag", 276288L, 53L)

# ---- Main extraction function ----------------------------------------------

#' Extract named parameters from CmdStan CSV files via column-pass.
#'
#' @param model_id       Character; used to discover CSV files if csv_files is NULL.
#' @param param_col_map  Named integer vector: param_name -> 0-indexed column position.
#' @param csv_files      Optional character vector of chain CSV paths (overrides discovery).
#' @param scratchpad     Directory for temporary per-chain output CSVs.
#' @return Data frame with columns .chain, .iteration, .draw, and one per param.
extract_stan_csv_params <- function(model_id = NULL,
                                    param_col_map,
                                    csv_files  = NULL,
                                    scratchpad = file.path(tempdir(), "stan_extract")) {
  # Locate Python script: src/extract_stan_csv_params.py relative to project root.
  # 'model_root' is set by bootstrap.R before sourcing this file.
  root <- tryCatch(model_root, error = function(e) getwd())
  py_script <- file.path(root, "src", "extract_stan_csv_params.py")
  if (!file.exists(py_script)) {
    py_script <- file.path(getwd(), "src", "extract_stan_csv_params.py")
  }
  if (!file.exists(py_script)) {
    stop("Cannot find src/extract_stan_csv_params.py. ",
         "Run from the project root or ensure 'paths$model' is set by bootstrap.R.")
  }

  if (is.null(csv_files)) {
    csv_files <- discover_cmdstan_csv_files(model_id)
    if (length(csv_files) == 0L) {
      stop("No CSV files found for model_id: ", model_id)
    }
  }

  dir.create(scratchpad, showWarnings = FALSE, recursive = TRUE)

  param_pairs <- paste0(names(param_col_map), ":", as.integer(param_col_map))

  chain_dfs <- lapply(seq_along(csv_files), function(chain_idx) {
    tmp_out <- file.path(scratchpad, sprintf("chain_%d_extract.csv", chain_idx))
    args    <- c(shQuote(py_script),
                 shQuote(csv_files[[chain_idx]]),
                 as.character(chain_idx),
                 shQuote(tmp_out),
                 param_pairs)
    ret <- system2("python3", args = args, stdout = FALSE, stderr = FALSE)
    if (ret != 0L || !file.exists(tmp_out)) {
      stop("Python extraction failed for chain ", chain_idx,
           " (CSV: ", csv_files[[chain_idx]], ")")
    }
    utils::read.csv(tmp_out, check.names = FALSE)
  })

  result  <- do.call(rbind, chain_dfs)
  n_iter  <- max(result$.iteration, na.rm = TRUE)
  result$.draw <- as.integer((result$.chain - 1L) * n_iter + result$.iteration)
  result
}
