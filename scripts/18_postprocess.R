# Stage 18 post-processing: nu summary, LOO comparison (Normal vs t-errors).
#
# Memory-efficient: extracts only log_lik columns from Stan CSV files using
# shell tools (grep + cut), never loads all 97K parameters into R.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/18_postprocess.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr", "loo", "data.table"))

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tibble); library(loo); library(data.table)
})

# ---- Helper: extract named columns from Stan CSV files via shell pipeline ---
# Never loads all parameters — uses grep+cut so only the requested columns
# enter R memory. Works for both log_lik extraction and scalar params like nu.

read_cols_from_csv <- function(csv_files, pattern) {
  header_fields <- strsplit(
    grep("^[^#]", readLines(csv_files[1], n = 200), value = TRUE)[1],
    ","
  )[[1]]
  idx <- grep(pattern, header_fields)
  if (length(idx) == 0) stop("No columns matching '", pattern, "' in ", csv_files[1])

  first_col <- min(idx)
  last_col  <- max(idx)
  message("  '", pattern, "' columns: ", first_col, "-", last_col,
          " (", length(idx), " cols)")

  chains <- lapply(csv_files, function(f) {
    cmd <- sprintf("grep -v '^#' %s | cut -d',' -f%d-%d",
                   shQuote(f), first_col, last_col)
    mat <- as.matrix(fread(cmd = cmd, header = TRUE))
    storage.mode(mat) <- "double"
    mat
  })
  do.call(rbind, chains)   # draws x cols
}

# ---- Posterior summary for nu (extracted from CSV, no full fit load) --------

message("Extracting nu draws from stage 18 CSVs...")
csv18   <- discover_cmdstan_csv_files("18_real_lowrank_a_diag_b_t")
nu_mat  <- read_cols_from_csv(csv18, "^nu$")
nu_vec  <- as.numeric(nu_mat)

nu_summary <- tibble::tibble(
  variable = "nu",
  mean     = mean(nu_vec),
  median   = median(nu_vec),
  sd       = sd(nu_vec),
  q5       = quantile(nu_vec, 0.05),
  q95      = quantile(nu_vec, 0.95)
)
cat("\nPosterior summary for nu (degrees of freedom):\n")
print(nu_summary)
rm(nu_mat, nu_vec); gc()

# ---- LOO for stage 18 -------------------------------------------------------

message("\nExtracting log_lik from stage 18 CSVs...")
ll_t  <- read_cols_from_csv(csv18, "^log_lik")
message("  Draws matrix: ", nrow(ll_t), " x ", ncol(ll_t))
loo_t <- loo(ll_t, cores = 1)
print(loo_t)
rm(ll_t); gc()

# ---- LOO for stage 10 -------------------------------------------------------

message("\nExtracting log_lik from stage 10 CSVs...")
csv10 <- discover_cmdstan_csv_files("10_real_lowrank_a_diag_b")
ll_n  <- read_cols_from_csv(csv10, "^log_lik")
message("  Draws matrix: ", nrow(ll_n), " x ", ncol(ll_n))
loo_n <- loo(ll_n, cores = 1)
print(loo_n)
rm(ll_n); gc()

# ---- Compare ----------------------------------------------------------------

comparison <- loo_compare(loo_n, loo_t)
cat("\nLOO comparison (Normal errors vs t-errors):\n")
print(comparison)

write.csv(as.data.frame(comparison),
          file.path(paths$tables, "18_loo_normal_vs_t.csv"))
message("Saved: ", file.path(paths$tables, "18_loo_normal_vs_t.csv"))

cat("\nnu posterior mean:  ", round(nu_summary$mean,   2), "\n")
cat("nu posterior median:", round(nu_summary$median, 2), "\n")
cat("nu 95% CI: [", round(nu_summary$q5, 2), ",", round(nu_summary$q95, 2), "]\n")

message("\nStage 18 post-processing complete.")
