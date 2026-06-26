# Stage 21: Fit the low-rank Sigma_a + diagonal Sigma_b model with K=4 factors
# and Student-t residuals.
#
# Extends stage 18 (K=2 t-errors) to K=4, motivated by rank selection results
# (stage 17): CV-ELPD is monotonically increasing through K=4 with no elbow.
#
# The t-errors model is used because it is clearly preferred over Gaussian
# (DELTA_ELPD = +11,668 at K=2; stage 18).
#
# Output:
#   fits/21_real_k4_t_fit.rds
#   outputs/tables/21_real_k4_t_posterior_summary.csv
#   outputs/tables/21_real_k4_t_loading_table.csv     (via write_fit_outputs)
#   outputs/figures/21_real_k4_t_chain_plots/

Sys.setenv(BFA_SIM_RANK_A  = "4",
           BFA_MAX_TREEDEPTH = "10",
           BFA_ADAPT_DELTA   = "0.92",
           BFA_ITER_WARMUP   = "1000",
           BFA_ITER_SAMPLING = "1000")

script_arg  <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "model/scripts/21_fit_k4_t.R"
source(file.path(dirname(dirname(normalizePath(script_path, mustWork = TRUE))), "src", "bootstrap.R"))
check_packages(c(required_base_packages, "cmdstanr"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

stopifnot(simulation$rank_a == 4L)
message("rank_a = ", simulation$rank_a)

model_objects_path  <- file.path(paths$processed, "model_objects.rds")
Y_path              <- file.path(paths$processed, "Y_scaled.csv")
stan_data_path_k4   <- file.path(paths$stan_data, "real_k4_t_stan_data.rds")

if (!file.exists(model_objects_path)) stop("Missing model_objects.rds. Run stage 02 first.")
if (!file.exists(Y_path))             stop("Missing Y_scaled.csv. Run stage 02 first.")

model_objects <- readRDS(model_objects_path)
Y_real        <- as.matrix(read_csv(Y_path, show_col_types = FALSE))

# Build K=4 stan data (same structure as K=2 but Q_a=4)
stan_data <- build_lowrank_a_diag_b_stan_data(
  model_objects = model_objects,
  Y             = Y_real,
  rank_a        = 4L
)
saveRDS(stan_data, stan_data_path_k4)
message("K=4 Stan data built: Q_a=", stan_data$Q_a,
        ", N_lambda_a_free=", stan_data$N_lambda_a_free)

# PCA warm-start (same logic as stage 18 but with K=4)
base_init_fn <- build_pca_init(stan_data)
init_fn <- function(chain_id = 1) {
  init <- base_init_fn(chain_id)
  init$mu <- NULL    # t-errors model has no mu
  init$nu <- 30.0    # start near Normal
  init
}

fit <- fit_cmdstan_model(
  stan_file  = file.path(paths$stan, "additive_lowrank_a_diag_b_t.stan"),
  stan_data  = stan_data,
  fit_path   = file.path(paths$fits, "21_real_k4_t_fit.rds"),
  seed       = pipeline$seed + 210L,
  model_id   = "21_real_k4_t",
  init       = init_fn
)

write_fit_outputs(fit, "21_real_k4_t")

if (!is.null(fit)) {
  nu_summary <- fit$summary("nu")
  cat("\nPosterior summary for nu (degrees of freedom):\n")
  print(nu_summary)

  write_notes(
    c(
      "Stage 21: K=4 factor model with Student-t residuals",
      "",
      paste0("N = ", nrow(Y_real), " player-season rows"),
      paste0("P = ", ncol(Y_real), " scaled features"),
      paste0("I = ", length(model_objects$player_levels), " unique players"),
      paste0("S = ", length(model_objects$season_levels), " seasons"),
      "Q_a = 4 player latent factors",
      "",
      "Student-t residuals: nu ~ Gamma(2, 0.1)",
      paste0("nu posterior mean: ", round(nu_summary$mean, 2)),
      paste0("nu 95% CI: [", round(nu_summary$q5, 2), ", ", round(nu_summary$q95, 2), "]"),
      "",
      "Motivation: stage 17 CV-ELPD shows monotonic improvement through K=4.",
      "Interpretation in stage 22."
    ),
    file.path(paths$notes, "21_real_k4_t_notes.txt")
  )
}

message("Stage 21 (K=4 t-errors fit) complete.")
