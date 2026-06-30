# Rank selection for the real player-season data.
#
# Each player-variable pair is assigned to one fold. When a fold is fitted,
# the code leaves every assigned variable out in every observed season for
# that player. Several variables for the same player may be left out together.

# -------------------------------------------------------------------------
# Prepare the data
# -------------------------------------------------------------------------

prepare_player_season_data = function(
    raw_data,
    player_column,
    season_column,
    variable_columns) {
  player_id = raw_data[[player_column]]
  season_id = raw_data[[season_column]]
  Y = as.matrix(raw_data[, variable_columns, drop = FALSE])
  storage.mode(Y) = "double"

  if (anyDuplicated(data.frame(player_id, season_id))) {
    stop("The data must have one row per observed player-season pair.")
  }
  if (any(!is.finite(Y))) {
    stop("The measurement columns must contain finite numeric values.")
  }

  center_y = colMeans(Y)
  scale_y = apply(Y, 2, sd)
  if (any(!is.finite(scale_y)) || any(scale_y <= 0)) {
    stop("Every measurement column must have a positive finite standard deviation.")
  }
  Y = sweep(sweep(Y, 2, center_y, "-"), 2, scale_y, "/")

  player_levels = unique(player_id)
  season_levels = unique(season_id)
  list(
    I = length(player_levels),
    J = length(season_levels),
    P = ncol(Y),
    N_obs = nrow(Y),
    player = match(player_id, player_levels),
    season = match(season_id, season_levels),
    Y = Y,
    player_levels = player_levels,
    season_levels = season_levels,
    variable_names = variable_columns,
    center_y = center_y,
    scale_y = scale_y
  )
}

# -------------------------------------------------------------------------
# Initial values
# -------------------------------------------------------------------------
# This adapts the PCA-based strategy in previous_lowrank_plus_diagonal/
# lowrank_diag.R. The only addition is that means are calculated after the
# current fold's held-out values have been removed.

column_sd = function(x, floor = 0.05) {
  out = apply(x, 2, sd, na.rm = TRUE)
  out[!is.finite(out)] = floor
  pmax(out, floor)
}

rough_additive_effects = function(dataset, split) {
  Y_train = matrix(NA_real_, dataset$N_obs, dataset$P)
  Y_train[split$train_index] = dataset$Y[split$train_index]

  # Rough season effects.
  B0 = matrix(0, dataset$J, dataset$P)
  for (season in seq_len(dataset$J)) {
    B0[season, ] = colMeans(
      Y_train[dataset$season == season, , drop = FALSE],
      na.rm = TRUE
    )
  }
  B0[!is.finite(B0)] = 0
  B0 = sweep(B0, 2, colMeans(B0), "-")

  # Rough player effects after removing the season effects.
  Y_without_season = Y_train - B0[dataset$season, ]
  A0 = matrix(0, dataset$I, dataset$P)
  for (player in seq_len(dataset$I)) {
    A0[player, ] = colMeans(
      Y_without_season[dataset$player == player, , drop = FALSE],
      na.rm = TRUE
    )
  }
  A0[!is.finite(A0)] = 0
  A0 = sweep(A0, 2, colMeans(A0), "-")

  error0 = Y_train - A0[dataset$player, ] - B0[dataset$season, ]
  list(
    A0 = A0,
    B0 = B0,
    sigma_b0 = column_sd(B0),
    sigma_e0 = column_sd(error0)
  )
}

make_diagonal_init = function(dataset, split) {
  start = rough_additive_effects(dataset, split)
  sigma_a0 = column_sd(start$A0)

  list(
    z_a = sweep(start$A0, 2, sigma_a0, "/"),
    z_b = sweep(start$B0, 2, start$sigma_b0, "/"),
    sigma_a = sigma_a0,
    sigma_b = start$sigma_b0,
    sigma_e = start$sigma_e0
  )
}

# Rotate the PCA loadings to the lower-triangular convention used in Stan.
make_plt = function(L) {
  K = ncol(L)
  decomposition = qr(t(L[seq_len(K), , drop = FALSE]))
  Q = qr.Q(decomposition)
  R = qr.R(decomposition)
  signs = sign(diag(R))
  signs[signs == 0] = 1
  L = L %*% Q %*% diag(signs, K)

  for (j in seq_len(K)) {
    for (h in seq_len(K)) {
      if (h > j) L[j, h] = 0
    }
  }
  for (h in seq_len(K)) {
    if (L[h, h] < 0) L[, h] = -L[, h]
  }
  L
}

make_lowrank_init = function(dataset, split, K) {
  start = rough_additive_effects(dataset, split)

  # PCA estimate of the loading matrix.
  player_covariance = crossprod(start$A0) / dataset$I
  decomposition = eigen(player_covariance, symmetric = TRUE)
  noise_variance = if (K < dataset$P) {
    mean(decomposition$values[(K + 1):dataset$P])
  } else {
    min(decomposition$values) * 0.5
  }
  loading_values = pmax(decomposition$values[seq_len(K)] - noise_variance, 1e-4)
  L0 = decomposition$vectors[, seq_len(K), drop = FALSE] %*%
    diag(sqrt(loading_values), K)
  L0 = make_plt(L0)

  # Starting factor scores and variable-specific player effects.
  psi_a0 = sqrt(pmax(diag(player_covariance) - rowSums(L0^2), 0.05^2))
  W = sweep(L0, 1, 1 / psi_a0^2, "*")
  F0 = start$A0 %*% W %*% solve(diag(K) + crossprod(L0, W))
  U0 = start$A0 - F0 %*% t(L0)
  psi_a0 = column_sd(U0)

  # Add a small perturbation, as in lowrank_diag.R, and store only the
  # lower-triangular entries sampled by the Stan model.
  L0 = L0 + matrix(rnorm(dataset$P * K, 0, 0.02), dataset$P, K)
  for (h in seq_len(K)) {
    for (j in seq_len(dataset$P)) {
      if (h > j) L0[j, h] = 0
    }
    L0[h, h] = abs(L0[h, h]) + 0.05
  }

  n_lower = K * dataset$P - (K * (K + 1)) %/% 2
  lambda_lower = numeric(n_lower)
  pos = 1
  for (h in seq_len(K)) {
    if (h < dataset$P) {
      for (j in (h + 1):dataset$P) {
        lambda_lower[pos] = L0[j, h]
        pos = pos + 1
      }
    }
  }

  list(
    F = F0,
    z_u = sweep(U0, 2, psi_a0, "/"),
    z_b = sweep(start$B0, 2, start$sigma_b0, "/"),
    psi_a = psi_a0,
    sigma_b = start$sigma_b0,
    sigma_e = start$sigma_e0,
    lambda_diag = diag(L0[seq_len(K), seq_len(K), drop = FALSE]),
    lambda_lower = lambda_lower
  )
}

# -------------------------------------------------------------------------
# Assign left-out player-variable pairs
# -------------------------------------------------------------------------

assign_folds = function(dataset, n_folds = 5, seed = 11) {
  if (length(n_folds) != 1 || n_folds < 2 || n_folds > dataset$I ||
      n_folds != as.integer(n_folds)) {
    stop("n_folds must be an integer between 2 and the number of players.")
  }

  set.seed(seed)
  folds = matrix(NA_integer_, dataset$I, dataset$P)

  # Each variable has its own random permutation of the players.
  for (variable in seq_len(dataset$P)) {
    folds[sample(seq_len(dataset$I)), variable] =
      rep(seq_len(n_folds), length.out = dataset$I)
  }
  folds
}

cell_index = function(mask) {
  out = which(mask, arr.ind = TRUE)
  storage.mode(out) = "integer"
  out
}

make_split = function(dataset, folds, fold) {
  # is_test has one row per observed player-season pair and one column per
  # variable. If (player, variable) is assigned to this fold, every season
  # observed for that player-variable pair is left out. Several variables for
  # the same player may be assigned to this fold.
  is_test = folds[dataset$player, , drop = FALSE] == fold
  list(
    train_index = cell_index(!is_test),
    test_index = cell_index(is_test)
  )
}

make_stan_data = function(dataset, split, rank) {
  out = list(
    N_obs = dataset$N_obs,
    I = dataset$I,
    J = dataset$J,
    P = dataset$P,
    player = dataset$player,
    season = dataset$season,
    N_train = nrow(split$train_index),
    train_obs = split$train_index[, 1],
    train_var = split$train_index[, 2],
    y_train = dataset$Y[split$train_index],
    N_test = nrow(split$test_index),
    test_obs = split$test_index[, 1],
    test_var = split$test_index[, 2],
    y_test = dataset$Y[split$test_index],
    scale_scale = 1
  )

  if (rank > 0) {
    out$K = rank
    out$lambda_scale = 1
  }
  out
}

# -------------------------------------------------------------------------
# Compile and fit the Stan models
# -------------------------------------------------------------------------

add_stan_runtime_to_path_on_windows = function() {
  if (.Platform$OS.type != "windows") {
    return(invisible(NULL))
  }

  # Stan executables need these C++ runtime libraries on Windows.
  candidates = c(
    file.path(
      cmdstanr::cmdstan_path(),
      "stan", "lib", "stan_math", "lib", "tbb"
    ),
    file.path(Sys.getenv("RTOOLS44_HOME", "C:/rtools44"), "ucrt64", "bin")
  )
  runtime_directories = candidates[dir.exists(candidates)]
  Sys.setenv(
    PATH = paste(c(runtime_directories, Sys.getenv("PATH")),
                 collapse = .Platform$path.sep)
  )
}

compile_models = function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Install cmdstanr and configure cmdstanr::set_cmdstan_path().")
  }
  add_stan_runtime_to_path_on_windows()

  list(
    diagonal = cmdstanr::cmdstan_model("stan/player_season_diagonal.stan"),
    lowrank = cmdstanr::cmdstan_model("stan/player_season_lowrank.stan")
  )
}

fit_one_rank = function(
    models,
    dataset,
    split,
    rank,
    seed,
    chains,
    parallel_chains,
    iter_warmup,
    iter_sampling,
    adapt_delta,
    refresh,
    method = "sample") {
  if (rank == 0) {
    model = models$diagonal
    initial_values = function(chain_id = 1) {
      make_diagonal_init(dataset, split)
    }
  } else {
    model = models$lowrank
    initial_values = function(chain_id = 1) {
      make_lowrank_init(dataset, split, rank)
    }
  }

  stan_data = make_stan_data(dataset, split, rank)

  if (method == "pathfinder") {
    # Pathfinder: fast approximate inference via normalizing flows.
    # Returns approximate posterior draws in minutes rather than hours.
    # Use for K-range screening only; confirm promising ranks with full NUTS.
    # num_paths=4 gives 4 independent runs (analogous to chains); draws per path
    # are aggregated. Note: Pathfinder quality degrades for very high-dimensional
    # posteriors; treat results as directional, not definitive.
    model$pathfinder(
      data       = stan_data,
      seed       = seed,
      num_paths  = max(chains, 4),
      draws      = iter_sampling,
      refresh    = refresh,
      init       = initial_values
    )
  } else {
    model$sample(
      data            = stan_data,
      seed            = seed,
      init            = initial_values,
      chains          = chains,
      parallel_chains = parallel_chains,
      iter_warmup     = iter_warmup,
      iter_sampling   = iter_sampling,
      adapt_delta     = adapt_delta,
      refresh         = refresh
    )
  }
}

# -------------------------------------------------------------------------
# Calculate held-out ELPD
# -------------------------------------------------------------------------

log_mean_exp = function(x) {
  maximum = max(x)
  maximum + log(mean(exp(x - maximum)))
}

score_fit = function(fit, split, rank, fold) {
  n_test = nrow(split$test_index)
  draw_names = paste0("log_lik_test[", seq_len(n_test), "]")
  log_lik = fit$draws(variables = "log_lik_test", format = "matrix")
  log_lik = log_lik[, draw_names, drop = FALSE]

  data.frame(
    rank = rank,
    fold = fold,
    log_score = apply(log_lik, 2, log_mean_exp)
  )
}

summarize_elpd = function(scores, ranks) {
  data.frame(
    rank = ranks,
    elpd = vapply(ranks, function(rank) {
      sum(scores$log_score[scores$rank == rank])
    }, numeric(1)),
    elpd_per_value = vapply(ranks, function(rank) {
      mean(scores$log_score[scores$rank == rank])
    }, numeric(1)),
    n_held_out_values = vapply(ranks, function(rank) {
      sum(scores$rank == rank)
    }, integer(1))
  )
}

# -------------------------------------------------------------------------
# Run the complete comparison
# -------------------------------------------------------------------------

run_rank_selection = function(
    dataset,
    ranks = 0:4,
    n_folds = 5,
    fold_seed = 11,
    fit_seed = 1000,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.98,
    refresh = 250,
    method = "sample") {
  if (any(ranks != as.integer(ranks)) || any(ranks < 0) ||
      any(ranks > dataset$P)) {
    stop("ranks must be integers between 0 and the number of variables.")
  }

  models = compile_models()
  folds = assign_folds(dataset, n_folds, fold_seed)
  scores = list()
  position = 1

  for (fold in seq_len(n_folds)) {
    split = make_split(dataset, folds, fold)

    for (rank in ranks) {
      message("Fitting fold ", fold, "/", n_folds, ", K = ", rank)
      fit = fit_one_rank(
        models, dataset, split, rank,
        seed = fit_seed + 100 * fold + rank,
        chains = chains,
        parallel_chains = parallel_chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        adapt_delta = adapt_delta,
        refresh = refresh,
        method = method
      )
      scores[[position]] = score_fit(fit, split, rank, fold)
      position = position + 1
    }
  }

  scores = do.call(rbind, scores)
  elpd_table = summarize_elpd(scores, ranks)
  holdouts = do.call(rbind, lapply(seq_len(n_folds), function(fold) {
    split = make_split(dataset, folds, fold)
    data.frame(
      fold = fold,
      held_out_player_variable_pairs = sum(folds == fold),
      held_out_values = nrow(split$test_index)
    )
  }))

  list(
    elpd = elpd_table,
    best_rank = elpd_table$rank[which.max(elpd_table$elpd)],
    holdouts = holdouts
  )
}
