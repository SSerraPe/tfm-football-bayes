# Simple Stan comparison of player-covariance models.
#
# K = 0 uses diagonal player, season, and error covariances.
# K > 0 uses a rank-K-plus-diagonal player covariance and diagonal season and
# error covariances.
#
# The main comparison assigns each player-variable pair to a fold. In a fitted
# fold, several variables for the same player may be omitted together across
# all observed seasons. The supplementary comparisons leave out a complete
# player-season row or every observation for a player.

`%||%` = function(x, y) {
  if (is.null(x)) y else x
}

# -------------------------------------------------------------------------
# Settings and simulated data
# -------------------------------------------------------------------------

default_sampling_config = function(
    chains = 4,
    parallel_chains = chains,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.95,
    max_treedepth = 12,
    refresh = 100) {
  list(
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = refresh
  )
}

make_unbalanced_design = function(I, J, min_seasons = 3, max_seasons = J) {
  repeat {
    design = do.call(
      rbind,
      lapply(seq_len(I), function(i) {
        data.frame(
          player = i,
          season = sort(sample(seq_len(J), sample(min_seasons:max_seasons, 1)))
        )
      })
    )

    if (all(tabulate(design$season, nbins = J) > 0)) {
      rownames(design) = NULL
      return(design)
    }
  }
}

make_true_loadings = function(P, K) {
  if (K == 0) {
    return(matrix(numeric(0), P, 0))
  }

  Lambda = matrix(0, P, K)
  for (h in seq_len(K)) {
    active = unique(c(h, seq(h, P, by = K)))
    Lambda[active, h] = rnorm(length(active), 0.75, 0.15)
    Lambda[h, h] = abs(Lambda[h, h]) + 0.35
  }

  if (K > 1) {
    shared = tail(seq_len(P), min(P, max(2, K + 1)))
    Lambda[shared, ] = Lambda[shared, ] +
      matrix(rnorm(length(shared) * K, 0.30, 0.08), length(shared), K)
  }

  Lambda
}

make_player_season_dataset = function(Y, player, season, standardize = TRUE) {
  Y = as.matrix(Y)
  storage.mode(Y) = "double"

  if (length(player) != nrow(Y) || length(season) != nrow(Y)) {
    stop("player and season must have one entry per row of Y.")
  }
  if (any(!is.finite(Y))) {
    stop("Y must contain finite values before fitting.")
  }

  if (standardize) {
    center_y = colMeans(Y)
    scale_y = apply(Y, 2, sd)
    if (any(!is.finite(scale_y)) || any(scale_y <= 0)) {
      stop("Every column of Y must have positive finite standard deviation.")
    }
    Y = sweep(sweep(Y, 2, center_y, "-"), 2, scale_y, "/")
  } else {
    center_y = rep(0, ncol(Y))
    scale_y = rep(1, ncol(Y))
  }

  player_levels = unique(player)
  season_levels = unique(season)

  list(
    I = length(player_levels),
    J = length(season_levels),
    P = ncol(Y),
    N_obs = nrow(Y),
    player = match(player, player_levels),
    season = match(season, season_levels),
    Y = Y,
    player_levels = player_levels,
    season_levels = season_levels,
    variable_names = colnames(Y),
    center_y = center_y,
    scale_y = scale_y
  )
}

simulate_player_season_data = function(
    I = 80,
    J = 6,
    P = 12,
    K_true = 2,
    min_seasons = 3,
    max_seasons = J,
    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  design = make_unbalanced_design(I, J, min_seasons, max_seasons)
  player = design$player
  season = design$season
  N_obs = nrow(design)

  Lambda = make_true_loadings(P, K_true)
  psi_a = runif(P, 0.25, 0.45)
  sigma_b = runif(P, 0.15, 0.30)
  sigma_e = runif(P, 0.40, 0.60)

  common_a = if (K_true == 0) {
    matrix(0, I, P)
  } else {
    matrix(rnorm(I * K_true), I, K_true) %*% t(Lambda)
  }

  U = sweep(matrix(rnorm(I * P), I, P), 2, psi_a, "*")
  B = sweep(matrix(rnorm(J * P), J, P), 2, sigma_b, "*")
  E = sweep(matrix(rnorm(N_obs * P), N_obs, P), 2, sigma_e, "*")
  Y_raw = common_a[player, ] + U[player, ] + B[season, ] + E

  dataset = make_player_season_dataset(Y_raw, player, season)
  Lambda_std = sweep(Lambda, 1, dataset$scale_y, "/")
  psi_a_std = psi_a / dataset$scale_y
  dataset$truth = list(
    K = K_true,
    Lambda = Lambda_std,
    psi_a = psi_a_std,
    Sigma_a = Lambda_std %*% t(Lambda_std) + diag(psi_a_std^2, P)
  )
  dataset
}

# -------------------------------------------------------------------------
# Leave-out schemes
# -------------------------------------------------------------------------

matrix_cell_index = function(mask) {
  index = which(mask, arr.ind = TRUE)
  if (length(index) == 0) {
    return(matrix(integer(0), 0, 2))
  }
  index
}

# Assign each player-variable combination to exactly one fold.
make_player_variable_folds = function(dataset, n_folds = 5, seed = 1) {
  set.seed(seed)
  folds = matrix(NA_integer_, dataset$I, dataset$P)

  for (p in seq_len(dataset$P)) {
    folds[sample(seq_len(dataset$I)), p] =
      rep(seq_len(n_folds), length.out = dataset$I)
  }

  folds
}

make_player_variable_split = function(dataset, folds, fold) {
  is_test = folds[dataset$player, , drop = FALSE] == fold
  list(
    train_index = matrix_cell_index(!is_test),
    test_index = matrix_cell_index(is_test)
  )
}

# Assign each observed player-season row to exactly one fold.
make_missing_season_folds = function(dataset, n_folds = 5, seed = 1) {
  set.seed(seed)
  folds = integer(dataset$N_obs)
  folds[sample(seq_len(dataset$N_obs))] =
    rep(seq_len(n_folds), length.out = dataset$N_obs)
  folds
}

make_missing_season_split = function(dataset, folds, fold) {
  is_test = matrix(folds == fold, dataset$N_obs, dataset$P)
  list(
    train_index = matrix_cell_index(!is_test),
    test_index = matrix_cell_index(is_test)
  )
}

# Assign each player to exactly one fold.
make_new_player_folds = function(dataset, n_folds = 5, seed = 1) {
  set.seed(seed)
  folds = integer(dataset$I)
  folds[sample(seq_len(dataset$I))] =
    rep(seq_len(n_folds), length.out = dataset$I)
  folds
}

make_new_player_split = function(dataset, folds, fold) {
  is_test = matrix(folds[dataset$player] == fold, dataset$N_obs, dataset$P)
  list(
    train_index = matrix_cell_index(!is_test),
    test_index = matrix_cell_index(is_test)
  )
}

stan_data_from_split = function(
    dataset,
    split,
    K = NULL,
    lambda_scale = 1,
    scale_scale = 1) {
  data = list(
    N_obs = dataset$N_obs,
    I = dataset$I,
    J = dataset$J,
    P = dataset$P,
    player = dataset$player,
    season = dataset$season,
    N_train = nrow(split$train_index),
    train_obs = as.integer(split$train_index[, 1]),
    train_var = as.integer(split$train_index[, 2]),
    y_train = as.numeric(dataset$Y[split$train_index]),
    N_test = nrow(split$test_index),
    test_obs = as.integer(split$test_index[, 1]),
    test_var = as.integer(split$test_index[, 2]),
    y_test = as.numeric(dataset$Y[split$test_index]),
    scale_scale = scale_scale
  )

  if (!is.null(K)) {
    data$K = K
    data$lambda_scale = lambda_scale
  }

  data
}

# -------------------------------------------------------------------------
# Stan compilation and initialization
# -------------------------------------------------------------------------

require_cmdstanr = function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Install cmdstanr and configure cmdstanr::set_cmdstan_path().")
  }
}

configure_cmdstan_windows_runtime = function() {
  if (.Platform$OS.type != "windows") {
    return(invisible(NULL))
  }

  tbb_dir = file.path(
    cmdstanr::cmdstan_path(),
    "stan", "lib", "stan_math", "lib", "tbb"
  )
  runtime_dirs = tbb_dir
  rtools_bin = file.path(Sys.getenv("RTOOLS44_HOME", "C:/rtools44"), "ucrt64", "bin")
  if (dir.exists(rtools_bin)) {
    runtime_dirs = c(runtime_dirs, rtools_bin)
  }

  if (any(!dir.exists(runtime_dirs))) {
    stop("The Stan runtime directory is missing. Check cmdstanr::cmdstan_path().")
  }

  Sys.setenv(
    PATH = paste(
      c(normalizePath(runtime_dirs, winslash = "/"), Sys.getenv("PATH")),
      collapse = .Platform$path.sep
    )
  )
}

compile_player_season_models = function(stan_dir = file.path("..", "stan")) {
  require_cmdstanr()
  configure_cmdstan_windows_runtime()

  list(
    diagonal = cmdstanr::cmdstan_model(
      file.path(stan_dir, "player_season_diagonal.stan")
    ),
    lowrank = cmdstanr::cmdstan_model(
      file.path(stan_dir, "player_season_lowrank.stan")
    )
  )
}

masked_training_matrix = function(dataset, split) {
  Y_train = matrix(NA_real_, dataset$N_obs, dataset$P)
  Y_train[split$train_index] = dataset$Y[split$train_index]
  Y_train
}

column_sd_with_floor = function(x, floor = 0.05) {
  out = apply(x, 2, sd, na.rm = TRUE)
  out[!is.finite(out)] = floor
  pmax(out, floor)
}

# Start from rough season and player means, following lowrank_diag.R.
make_additive_start = function(dataset, split, scale_floor = 0.05) {
  Y_train = masked_training_matrix(dataset, split)

  B0 = matrix(0, dataset$J, dataset$P)
  for (j in seq_len(dataset$J)) {
    B0[j, ] = colMeans(
      Y_train[dataset$season == j, , drop = FALSE],
      na.rm = TRUE
    )
  }
  B0[!is.finite(B0)] = 0
  B0 = sweep(B0, 2, colMeans(B0), "-")

  Y_without_season = Y_train - B0[dataset$season, ]
  A0 = matrix(0, dataset$I, dataset$P)
  for (i in seq_len(dataset$I)) {
    A0[i, ] = colMeans(
      Y_without_season[dataset$player == i, , drop = FALSE],
      na.rm = TRUE
    )
  }
  A0[!is.finite(A0)] = 0
  A0 = sweep(A0, 2, colMeans(A0), "-")

  E0 = Y_train - A0[dataset$player, ] - B0[dataset$season, ]

  list(
    A0 = A0,
    B0 = B0,
    sigma_b0 = column_sd_with_floor(B0, scale_floor),
    sigma_e0 = column_sd_with_floor(E0, scale_floor)
  )
}

# Rotate PCA loadings to the lower-triangular convention used by the Stan model.
make_plt = function(L) {
  K = ncol(L)
  qr_A = qr(t(L[seq_len(K), , drop = FALSE]))
  Q = qr.Q(qr_A)
  R = qr.R(qr_A)
  signs = sign(diag(R))
  signs[signs == 0] = 1
  L = L %*% Q %*% diag(signs, K)

  for (row in seq_len(K)) {
    for (column in seq_len(K)) {
      if (column > row) {
        L[row, column] = 0
      }
    }
  }

  for (column in seq_len(K)) {
    if (L[column, column] < 0) {
      L[, column] = -L[, column]
    }
  }

  L
}

make_diagonal_init = function(dataset, split, scale_floor = 0.05) {
  start = make_additive_start(dataset, split, scale_floor)
  sigma_a0 = column_sd_with_floor(start$A0, scale_floor)

  list(
    z_a = sweep(start$A0, 2, sigma_a0, "/"),
    z_b = sweep(start$B0, 2, start$sigma_b0, "/"),
    sigma_a = sigma_a0,
    sigma_b = start$sigma_b0,
    sigma_e = start$sigma_e0
  )
}

make_lowrank_init = function(dataset, split, K, scale_floor = 0.05) {
  start = make_additive_start(dataset, split, scale_floor)
  S_A = crossprod(start$A0) / dataset$I
  eig = eigen(S_A, symmetric = TRUE)
  sigma2_iso = if (K < dataset$P) {
    mean(eig$values[(K + 1):dataset$P])
  } else {
    min(eig$values) * 0.5
  }

  load_values = pmax(eig$values[seq_len(K)] - sigma2_iso, 1e-4)
  L0 = eig$vectors[, seq_len(K), drop = FALSE] %*%
    diag(sqrt(load_values), K)
  L0 = make_plt(L0)

  psi_a0 = sqrt(pmax(diag(S_A) - rowSums(L0^2), scale_floor^2))
  W = sweep(L0, 1, 1 / psi_a0^2, "*")
  F0 = start$A0 %*% W %*% solve(diag(K) + crossprod(L0, W))
  U0 = start$A0 - F0 %*% t(L0)
  psi_a0 = column_sd_with_floor(U0, scale_floor)

  lambda_lower = numeric(K * dataset$P - (K * (K + 1)) %/% 2)
  pos = 1
  for (column in seq_len(K)) {
    if (column < dataset$P) {
      for (row in (column + 1):dataset$P) {
        lambda_lower[pos] = L0[row, column]
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

sample_rank_model = function(
    compiled_models,
    dataset,
    split,
    rank,
    seed,
    sampling_config = default_sampling_config()) {
  if (rank == 0) {
    model = compiled_models$diagonal
    data = stan_data_from_split(dataset, split)
    init = function(chain_id = 1) make_diagonal_init(dataset, split)
  } else {
    model = compiled_models$lowrank
    data = stan_data_from_split(dataset, split, K = rank)
    init = function(chain_id = 1) make_lowrank_init(dataset, split, rank)
  }

  fit = do.call(
    model$sample,
    c(list(data = data, seed = seed, init = init), sampling_config)
  )

  if (any(fit$return_codes() != 0)) {
    stop("Stan chain failure:\n", paste(fit$output(), collapse = "\n"))
  }

  fit
}

# -------------------------------------------------------------------------
# Held-out ELPD
# -------------------------------------------------------------------------

draw_matrix = function(fit, variable, length) {
  variables = paste0(variable, "[", seq_len(length), "]")
  draws = fit$draws(variables = variable, format = "matrix")
  draws[, variables, drop = FALSE]
}

log_mean_exp_columns = function(x) {
  apply(x, 2, function(values) {
    maximum = max(values)
    maximum + log(mean(exp(values - maximum)))
  })
}

extract_cell_scores = function(fit, split, rank, fold) {
  data.frame(
    rank = rank,
    fold = fold,
    log_score = log_mean_exp_columns(
      draw_matrix(fit, "log_lik_test", nrow(split$test_index))
    )
  )
}

summarize_elpd_scores = function(scores) {
  do.call(
    rbind,
    lapply(split(scores, scores$rank), function(x) {
      data.frame(
        rank = x$rank[1],
        elpd = sum(x$log_score),
        elpd_per_value = mean(x$log_score),
        n_held_out_values = nrow(x)
      )
    })
  )
}

run_bayesian_rank_cv = function(
    dataset,
    fit_ranks = 0:4,
    n_folds = 5,
    fold_seed = 11,
    fit_seed = 1000,
    compiled_models = NULL,
    sampling_config = default_sampling_config()) {
  models = compiled_models %||% compile_player_season_models()
  folds = make_player_variable_folds(dataset, n_folds, fold_seed)
  scores = list()
  pos = 1

  for (fold in seq_len(n_folds)) {
    split = make_player_variable_split(dataset, folds, fold)

    for (rank in fit_ranks) {
      message("Fitting player-variable fold ", fold, "/", n_folds,
              ", K = ", rank)
      fit = sample_rank_model(
        models, dataset, split, rank,
        seed = fit_seed + 100 * fold + rank,
        sampling_config = sampling_config
      )
      scores[[pos]] = extract_cell_scores(fit, split, rank, fold)
      pos = pos + 1
    }
  }

  scores = do.call(rbind, scores)
  holdouts = do.call(
    rbind,
    lapply(seq_len(n_folds), function(fold) {
      split = make_player_variable_split(dataset, folds, fold)
      data.frame(
        fold = fold,
        held_out_combinations = sum(folds == fold),
        held_out_values = nrow(split$test_index)
      )
    })
  )

  summary = summarize_elpd_scores(scores)
  list(
    scores = scores,
    summary = summary,
    best_rank = summary$rank[which.max(summary$elpd)],
    holdouts = holdouts
  )
}

# -------------------------------------------------------------------------
# Supplementary leave-out comparisons
# -------------------------------------------------------------------------

extract_grouped_scores = function(fit, dataset, split, rank, fold, groups) {
  test_index = split$test_index
  log_lik = draw_matrix(fit, "log_lik_test", nrow(test_index))
  group_for_value = as.character(groups[test_index[, 1]])

  do.call(
    rbind,
    lapply(unique(group_for_value), function(group) {
      values = which(group_for_value == group)
      joint_log_lik = rowSums(log_lik[, values, drop = FALSE])
      maximum = max(joint_log_lik)

      data.frame(
        rank = rank,
        fold = fold,
        group = group,
        n_values = length(values),
        log_score = maximum + log(mean(exp(joint_log_lik - maximum)))
      )
    })
  )
}

summarize_grouped_scores = function(scores) {
  do.call(
    rbind,
    lapply(split(scores, scores$rank), function(x) {
      data.frame(
        rank = x$rank[1],
        elpd = sum(x$log_score),
        n_left_out_units = nrow(x),
        n_held_out_values = sum(x$n_values)
      )
    })
  )
}

run_bayesian_supplementary_cv = function(
    dataset,
    task = c("missing_season", "new_player"),
    fit_ranks = 0:4,
    n_folds = 5,
    fold_seed = 21,
    fit_seed = 3000,
    compiled_models = NULL,
    sampling_config = default_sampling_config()) {
  task = match.arg(task)
  models = compiled_models %||% compile_player_season_models()

  if (task == "missing_season") {
    folds = make_missing_season_folds(dataset, n_folds, fold_seed)
    split_function = make_missing_season_split
    groups = seq_len(dataset$N_obs)
  } else {
    folds = make_new_player_folds(dataset, n_folds, fold_seed)
    split_function = make_new_player_split
    groups = dataset$player
  }

  scores = list()
  pos = 1
  for (fold in seq_len(n_folds)) {
    split = split_function(dataset, folds, fold)

    for (rank in fit_ranks) {
      message("Fitting ", task, " fold ", fold, "/", n_folds,
              ", K = ", rank)
      fit = sample_rank_model(
        models, dataset, split, rank,
        seed = fit_seed + 100 * fold + rank,
        sampling_config = sampling_config
      )
      scores[[pos]] =
        extract_grouped_scores(fit, dataset, split, rank, fold, groups)
      pos = pos + 1
    }
  }

  scores = do.call(rbind, scores)
  summary = summarize_grouped_scores(scores)
  list(
    scores = scores,
    summary = summary,
    best_rank = summary$rank[which.max(summary$elpd)]
  )
}

# -------------------------------------------------------------------------
# Simulation study used in the note
# -------------------------------------------------------------------------

run_bayesian_simulation_study = function(
    scenarios = 0:3,
    replicates = 2,
    fit_ranks = 0:4,
    I = 80,
    J = 6,
    P = 12,
    n_folds = 2,
    seed = 500,
    sampling_config = default_sampling_config()) {
  models = compile_player_season_models()
  runs = list()
  summaries = list()
  holdouts = list()
  pos = 1

  for (true_rank in scenarios) {
    for (replicate in seq_len(replicates)) {
      message("Simulation truth K = ", true_rank, ", replicate = ", replicate)
      dataset = simulate_player_season_data(
        I, J, P, true_rank,
        seed = seed + 10000 * true_rank + replicate
      )
      cv = run_bayesian_rank_cv(
        dataset = dataset,
        fit_ranks = fit_ranks,
        n_folds = n_folds,
        fold_seed = seed + replicate,
        fit_seed = seed + 1000 * replicate,
        compiled_models = models,
        sampling_config = sampling_config
      )

      runs[[pos]] = data.frame(
        true_rank = true_rank,
        replicate = replicate,
        best_rank = cv$best_rank
      )
      summaries[[pos]] = transform(
        cv$summary,
        true_rank = true_rank,
        replicate = replicate
      )
      holdouts[[pos]] = transform(
        cv$holdouts,
        true_rank = true_rank,
        replicate = replicate
      )
      pos = pos + 1
    }
  }

  list(
    runs = do.call(rbind, runs),
    cv_summary = do.call(rbind, summaries),
    holdouts = do.call(rbind, holdouts)
  )
}
