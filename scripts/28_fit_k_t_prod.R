# Stage 28: Fit the t-errors model at a chosen rank K (from BFA_SIM_RANK_A).
#
# Use this when Task O's Pathfinder screening selects K != 2.  Builds a fresh
# stan_data at the chosen K and fits additive_lowrank_a_diag_b_t.stan.
#
# Usage:
#   BFA_RUN_STAN=true BFA_REUSE_FIT=false BFA_SIM_RANK_A=<K> Rscript scripts/28_fit_k_t_prod.R
#
# If BFA_SIM_RANK_A=2, the result is identical to stage 18; use stage 18 instead.
# This script is the Task O production runner for K > 2.

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/28_fit_k_t_prod.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

rank_a   <- simulation$rank_a   # reads BFA_SIM_RANK_A; default 2
model_id <- paste0("28_real_lowrank_a_diag_b_t_k", rank_a)
message("Task O fit: t-errors model at K=", rank_a, "  (model_id: ", model_id, ")")

# ---- Load rebuilt model objects and Y ----------------------------------------

model_objects_path <- file.path(paths$processed, "model_objects.rds")
if (!file.exists(model_objects_path))
  stop("Missing model_objects.rds. Run stage 02 first.", call. = FALSE)

model_objects <- readRDS(model_objects_path)
Y             <- model_objects$Y

message("Data: N=", nrow(Y), " P=", ncol(Y), " I=",
        length(model_objects$player_levels), " S=", length(model_objects$season_levels))

# ---- Build stan_data at rank K -----------------------------------------------

stan_data <- build_lowrank_a_diag_b_stan_data(model_objects, Y, rank_a)
stan_data$compute_log_lik <- 0L   # skip log_lik loop on the fitting run

saveRDS(stan_data,
        file.path(paths$stan_data,
                  paste0("real_lowrank_a_diag_b_k", rank_a, "_stan_data.rds")))

# ---- PCA init adapted for t-model -------------------------------------------

base_init_fn <- build_pca_init(stan_data)

init_fn <- function(chain_id = 1) {
  init     <- base_init_fn(chain_id)
  init$mu  <- NULL   # t-model has no mu
  init$nu  <- 30.0   # start near Normal
  init
}

# ---- Fit ---------------------------------------------------------------------

fit_path <- file.path(paths$fits, paste0(model_id, "_fit.rds"))

fit <- fit_cmdstan_model(
  stan_file  = file.path(paths$stan, "additive_lowrank_a_diag_b_t.stan"),
  stan_data  = stan_data,
  fit_path   = fit_path,
  seed       = pipeline$seed + 280L,
  model_id   = model_id,
  init       = init_fn
)

write_fit_outputs(fit, model_id)

# ---- nu summary --------------------------------------------------------------

if (!is.null(fit)) {
  nu_sum <- tryCatch(fit$summary("nu"), error = function(e) NULL)
  if (!is.null(nu_sum)) {
    message(sprintf("\nnu posterior: mean=%.3f  95%%CI=[%.3f, %.3f]",
                    nu_sum$mean, nu_sum$q5, nu_sum$q95))
  }
  message("\nCSV files: ", paste(fit$output_files(), collapse = "\n  "))
}

message("\nStage 28 complete (K=", rank_a, ").")
message("To compute log_lik for LOO, call fit$generate_quantities() with ",
        "compute_log_lik=1 on the stan_data.")
