# Helpers for low-rank-plus-diagonal additive simulation.

center_columns <- function(x) {
  sweep(x, 2, colMeans(x), "-")
}

make_anchor_loadings <- function(P, Q, diag_min, diag_max, free_sd) {
  if (Q > P) {
    stop("Simulation rank cannot exceed the number of variables.", call. = FALSE)
  }
  L <- matrix(0, nrow = P, ncol = Q)
  diag_values <- stats::runif(Q, diag_min, diag_max)
  for (q in seq_len(Q)) {
    L[q, q] <- diag_values[q]
  }
  if (P > Q) {
    L[(Q + 1L):P, ] <- matrix(stats::rnorm((P - Q) * Q, 0, free_sd), nrow = P - Q)
  }
  L
}

anchor_free_values <- function(L, Q) {
  P <- nrow(L)
  if (P <= Q) return(numeric(0))
  as.vector(L[(Q + 1L):P, , drop = FALSE])
}

matrix_to_long_table <- function(x, matrix_name, row_names = NULL, col_prefix = "factor") {
  if (is.null(row_names)) row_names <- as.character(seq_len(nrow(x)))
  out <- expand.grid(
    observed_index = seq_len(nrow(x)),
    factor_index = seq_len(ncol(x)),
    KEEP.OUT.ATTRS = FALSE
  )
  out$value <- as.numeric(x[cbind(out$observed_index, out$factor_index)])
  tibble::as_tibble(out) |>
    dplyr::mutate(
      parameter = matrix_name,
      observed_variable = row_names[observed_index],
      factor = paste0(col_prefix, factor_index)
    ) |>
    dplyr::select(parameter, observed_index, observed_variable, factor_index, factor, value)
}

simulate_lowrank_additive <- function(model_objects, rank_a, rank_b, seed) {
  if (rank_a <= rank_b) {
    stop("Simulation requires rank_a > rank_b.", call. = FALSE)
  }

  set.seed(seed)

  Y_template <- model_objects$Y
  N <- nrow(Y_template)
  P <- ncol(Y_template)
  I <- length(model_objects$player_levels)
  S <- length(model_objects$season_levels)

  Lambda_a <- make_anchor_loadings(P, rank_a, diag_min = 0.45, diag_max = 0.85, free_sd = 0.20)
  Lambda_b <- make_anchor_loadings(P, rank_b, diag_min = 0.25, diag_max = 0.55, free_sd = 0.12)

  sigma_unique_a <- stats::runif(P, 0.15, 0.35)
  sigma_unique_b <- stats::runif(P, 0.08, 0.25)
  sigma_epsilon <- stats::runif(P, 0.35, 0.65)
  mu <- stats::rnorm(P, 0, 0.15)

  eta_a <- center_columns(matrix(stats::rnorm(I * rank_a), nrow = I, ncol = rank_a))
  eta_b <- center_columns(matrix(stats::rnorm(S * rank_b), nrow = S, ncol = rank_b))
  unique_a <- center_columns(matrix(stats::rnorm(I * P), nrow = I, ncol = P))
  unique_b <- center_columns(matrix(stats::rnorm(S * P), nrow = S, ncol = P))

  a <- eta_a %*% t(Lambda_a) + sweep(unique_a, 2, sigma_unique_a, `*`)
  b <- eta_b %*% t(Lambda_b) + sweep(unique_b, 2, sigma_unique_b, `*`)
  a <- center_columns(a)
  b <- center_columns(b)

  Y_sim <- matrix(NA_real_, nrow = N, ncol = P)
  for (n in seq_len(N)) {
    mean_n <- mu + a[model_objects$player_index[n], ] + b[model_objects$season_index[n], ]
    Y_sim[n, ] <- mean_n + stats::rnorm(P, 0, sigma_epsilon)
  }
  colnames(Y_sim) <- model_objects$variable_names

  Sigma_a <- Lambda_a %*% t(Lambda_a) + diag(sigma_unique_a^2, P)
  Sigma_b <- Lambda_b %*% t(Lambda_b) + diag(sigma_unique_b^2, P)
  Sigma_epsilon <- diag(sigma_epsilon^2, P)

  truth <- list(
    seed = seed,
    rank_a = rank_a,
    rank_b = rank_b,
    mu = mu,
    Lambda_a = Lambda_a,
    Lambda_b = Lambda_b,
    Psi_a = sigma_unique_a^2,
    Psi_b = sigma_unique_b^2,
    sigma_unique_a = sigma_unique_a,
    sigma_unique_b = sigma_unique_b,
    sigma_epsilon = sqrt(diag(Sigma_epsilon)),
    Sigma_epsilon = Sigma_epsilon,
    Sigma_a = Sigma_a,
    Sigma_b = Sigma_b,
    Sigma_a_diag = diag(Sigma_a),
    Sigma_b_diag = diag(Sigma_b),
    a = a,
    b = b,
    eta_a = eta_a,
    eta_b = eta_b
  )

  list(Y = Y_sim, truth = truth)
}
