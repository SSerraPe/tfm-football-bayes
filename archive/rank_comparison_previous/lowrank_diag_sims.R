library(cmdstanr)

set.seed(1)

# -----------------------------
# Settings
# -----------------------------

I = 100
J = 8
P = 20
K_true = 2
K_values = 1:4

iter_warmup = 500
iter_sampling = 500
chains = 4
n_check = 600

# -----------------------------
# Unbalanced player-season design
# -----------------------------

make_design = function(I, J, min_seasons = 3, max_seasons = 8) {
  out = vector("list", I)
  
  for (i in 1:I) {
    m_i = sample(min_seasons:max_seasons, 1)
    seasons_i = sort(sample(1:J, m_i))
    out[[i]] = data.frame(player = i, season = seasons_i)
  }
  
  do.call(rbind, out)
}

repeat {
  design = make_design(I, J, min_seasons = 3, max_seasons = J)
  if (all(tabulate(design$season, nbins = J)>0)) break
}

player = design$player
season = design$season
N_obs = nrow(design)

cat("Number of observations:", N_obs, "\n")
cat("Seasons per player:\n")
print(summary(tabulate(player, nbins = I)))
cat("Players per season:\n")
print(tabulate(season, nbins = J))

# -----------------------------
# Simulate crossed factor model
# -----------------------------

Lambda_true = matrix(0, P, K_true)

Lambda_true[1, 1] = 1.2
Lambda_true[2, 1] = 0.4
Lambda_true[2, 2] = 1.0

active1 = 3:8
active2 = 9:14
both = 15:18

Lambda_true[active1, 1] = rnorm(length(active1), 0.8, 0.2)
Lambda_true[active2, 2] = rnorm(length(active2), -0.8, 0.2)
Lambda_true[both, 1] = rnorm(length(both), 0.4, 0.1)
Lambda_true[both, 2] = rnorm(length(both), 0.4, 0.1)

psi_a_true = runif(P, 0.25, 0.45)
sigma_b_true = runif(P, 0.15, 0.30)
sigma_e_true = runif(P, 0.40, 0.60)

F_true = matrix(rnorm(I*K_true), I, K_true)

U_true = matrix(rnorm(I*P), I, P)
U_true = sweep(U_true, 2, psi_a_true, "*")

B_season_true = matrix(rnorm(J*P), J, P)
B_season_true = sweep(B_season_true, 2, sigma_b_true, "*")

E_true = matrix(rnorm(N_obs*P), N_obs, P)
E_true = sweep(E_true, 2, sigma_e_true, "*")

A_true = F_true%*%t(Lambda_true)+U_true
Y_raw = A_true[player, ]+B_season_true[season, ]+E_true

# -----------------------------
# Standardize observed data
# -----------------------------

center_y = colMeans(Y_raw)
scale_y = apply(Y_raw, 2, sd)

Y = sweep(Y_raw, 2, center_y, "-")
Y = sweep(Y, 2, scale_y, "/")

Lambda_true_std = sweep(Lambda_true, 1, scale_y, "/")
psi_a_true_std = psi_a_true/scale_y
sigma_b_true_std = sigma_b_true/scale_y
sigma_e_true_std = sigma_e_true/scale_y

B_a_true = Lambda_true_std%*%t(Lambda_true_std)
Sigma_a_true = B_a_true+diag(psi_a_true_std^2, P)
R_a_true = cov2cor(Sigma_a_true)

Psi_a_true = psi_a_true_std^2
Db_true = sigma_b_true_std^2
De_true = sigma_e_true_std^2
pop_var_true = diag(Sigma_a_true)+Db_true+De_true

# -----------------------------
# Stan model
# -----------------------------

stan_code = "
data {
  int<lower=1> N_obs;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> P;
  int<lower=1, upper=P> K;

  array[N_obs] int<lower=1, upper=I> player;
  array[N_obs] int<lower=1, upper=J> season;

  matrix[N_obs, P] Y;

  real<lower=0> lambda_scale;
  real<lower=0> scale_scale;
}

transformed data {
  int<lower=0> n_lower = K * P - (K * (K + 1)) %/% 2;
}

parameters {
  matrix[I, K] F;
  matrix[I, P] z_u;
  matrix[J, P] z_b;

  vector<lower=1e-6>[P] psi_a;
  vector<lower=1e-6>[P] sigma_b;
  vector<lower=1e-6>[P] sigma_e;

  vector<lower=0>[K] lambda_diag;
  vector[n_lower] lambda_lower;
}

transformed parameters {
  matrix[P, K] Lambda;
  vector<lower=0>[K] lambda_col_norm;

  Lambda = rep_matrix(0, P, K);

  {
    int pos;
    pos = 1;

    for (h in 1:K) {
      Lambda[h, h] = lambda_diag[h];

      if (h < P) {
        for (j in (h + 1):P) {
          Lambda[j, h] = lambda_lower[pos];
          pos += 1;
        }
      }
    }
  }

  for (h in 1:K) {
    lambda_col_norm[h] = sqrt(dot_self(Lambda[, h]));
  }
}

model {
  matrix[I, P] A_low;
  vector[N_obs] mu;

  to_vector(F) ~ std_normal();
  to_vector(z_u) ~ std_normal();
  to_vector(z_b) ~ std_normal();

  lambda_diag ~ normal(0, lambda_scale);
  lambda_lower ~ normal(0, lambda_scale);

  psi_a ~ student_t(3, 0, scale_scale);
  sigma_b ~ student_t(3, 0, scale_scale);
  sigma_e ~ student_t(3, 0, scale_scale);

  A_low = F * Lambda';

  for (p in 1:P) {
    for (n in 1:N_obs) {
      mu[n] =
        A_low[player[n], p]
        + psi_a[p] * z_u[player[n], p]
        + sigma_b[p] * z_b[season[n], p];
    }

    Y[, p] ~ normal(mu, sigma_e[p]);
  }
}
"

stan_file = write_stan_file(stan_code)

mod = cmdstan_model(stan_file, compile = FALSE)
mod$check_syntax()
mod$compile()

# -----------------------------
# Initialization
# -----------------------------

make_plt = function(L) {
  K = ncol(L)
  A = L[1:K,, drop = FALSE]
  
  qr_A = qr(t(A))
  Q = qr.Q(qr_A)
  R = qr.R(qr_A)
  
  sgn = sign(diag(R))
  sgn[sgn==0] = 1
  
  Q = Q%*%diag(sgn, K)
  L_tri = L%*%Q
  
  for (j in 1:K) {
    for (h in 1:K) {
      if (h>j) L_tri[j, h] = 0
    }
  }
  
  for (h in 1:K) {
    if (L_tri[h, h]<0) L_tri[, h] = -L_tri[, h]
  }
  
  L_tri
}

make_crossed_init = function(Y, player, season, I, J, K,
                             lambda_jitter = 0.02,
                             scale_floor = 0.05) {
  P = ncol(Y)
  
  season_counts = tabulate(season, nbins = J)
  B0 = sweep(rowsum(Y, season), 1, season_counts, "/")
  B0 = sweep(B0, 2, colMeans(B0), "-")
  
  Y_adj = Y-B0[season, ]
  
  player_counts = tabulate(player, nbins = I)
  A0 = sweep(rowsum(Y_adj, player), 1, player_counts, "/")
  A0 = sweep(A0, 2, colMeans(A0), "-")
  
  S_A = crossprod(A0)/I
  eig = eigen(S_A, symmetric = TRUE)
  
  vals = eig$values
  vecs = eig$vectors
  
  if (K<P) {
    sigma2_iso = mean(vals[(K+1):P])
  } else {
    sigma2_iso = min(vals)*0.5
  }
  
  load_vals = pmax(vals[1:K]-sigma2_iso, 1e-4)
  L0 = vecs[, 1:K, drop = FALSE]%*%diag(sqrt(load_vals), K)
  L0 = make_plt(L0)
  
  resid_var0 = pmax(diag(S_A)-rowSums(L0^2), scale_floor^2)
  psi_a0 = sqrt(resid_var0)
  
  psi_inv = 1/psi_a0^2
  W = sweep(L0, 1, psi_inv, "*")
  V = solve(diag(K)+crossprod(L0, W))
  F0 = A0%*%W%*%V
  
  U0 = A0-F0%*%t(L0)
  
  psi_a0 = pmax(apply(U0, 2, sd), scale_floor)
  z_u0 = sweep(U0, 2, psi_a0, "/")
  
  sigma_b0 = pmax(apply(B0, 2, sd), scale_floor)
  z_b0 = sweep(B0, 2, sigma_b0, "/")
  
  Ahat0 = F0%*%t(L0)+U0
  E0 = Y-Ahat0[player, ]-B0[season, ]
  sigma_e0 = pmax(apply(E0, 2, sd), scale_floor)
  
  L0 = L0+matrix(rnorm(P*K, 0, lambda_jitter), P, K)
  
  for (h in 1:K) {
    for (j in 1:P) {
      if (h>j) L0[j, h] = 0
    }
    L0[h, h] = abs(L0[h, h])+0.05
  }
  
  n_lower = K*P-(K*(K+1))%/%2
  
  lambda_diag = diag(L0[1:K, 1:K, drop = FALSE])
  lambda_lower = numeric(n_lower)
  
  pos = 1
  for (h in 1:K) {
    if (h<P) {
      for (j in (h+1):P) {
        lambda_lower[pos] = L0[j, h]
        pos = pos+1
      }
    }
  }
  
  list(
    F = F0,
    z_u = z_u0,
    z_b = z_b0,
    psi_a = psi_a0,
    sigma_b = sigma_b0,
    sigma_e = sigma_e0,
    lambda_diag = lambda_diag,
    lambda_lower = lambda_lower
  )
}

# -----------------------------
# Fit and check one rank
# -----------------------------

fit_one_rank = function(K_fit) {
  
  cat("\n\n=============================\n")
  cat("Fitting K =", K_fit, "\n")
  cat("=============================\n\n")
  
  data_list = list(
    N_obs = N_obs,
    I = I,
    J = J,
    P = P,
    K = K_fit,
    player = player,
    season = season,
    Y = Y,
    lambda_scale = 1,
    scale_scale = 1
  )
  
  make_init = function(chain_id = 1) {
    make_crossed_init(Y, player, season, I, J, K_fit)
  }
  
  time_fit = system.time({
    fit = mod$sample(
      data = data_list,
      seed = 1000+K_fit,
      chains = chains,
      parallel_chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = 0.95,
      max_treedepth = 12,
      init = make_init,
      refresh = 500
    )
  })
  
  print(time_fit)
  print(fit$diagnostic_summary())
  
  diag_summ = fit$summary(c("lp__", "lambda_col_norm", "psi_a", "sigma_b", "sigma_e"))
  
  lambda_arr = fit$draws(variables = "Lambda", format = "draws_array")
  psi_arr = fit$draws(variables = "psi_a", format = "draws_array")
  sigb_arr = fit$draws(variables = "sigma_b", format = "draws_array")
  sige_arr = fit$draws(variables = "sigma_e", format = "draws_array")
  lp_arr = fit$draws(variables = "lp__", format = "draws_array")
  
  n_iter = dim(lambda_arr)[1]
  n_chains = dim(lambda_arr)[2]
  all_draws = expand.grid(iter = 1:n_iter, chain = 1:n_chains)
  
  set.seed(2000+K_fit)
  use_draws = all_draws[sample(seq_len(nrow(all_draws)), min(n_check, nrow(all_draws))), ]
  n_use = nrow(use_draws)
  
  Lambda_vars = as.vector(outer(1:P, 1:K_fit, function(j, h) paste0("Lambda[", j, ",", h, "]")))
  psi_vars = paste0("psi_a[", 1:P, "]")
  sigb_vars = paste0("sigma_b[", 1:P, "]")
  sige_vars = paste0("sigma_e[", 1:P, "]")
  
  Lambda_mean = matrix(0, P, K_fit)
  B_a_mean = matrix(0, P, P)
  Sigma_a_mean = matrix(0, P, P)
  R_a_mean = matrix(0, P, P)
  
  psi2_mean = numeric(P)
  sigb2_mean = numeric(P)
  sige2_mean = numeric(P)
  pop_var_mean = numeric(P)
  
  eig_common = matrix(NA, n_use, K_fit)
  
  Q_true = qr.Q(qr(Lambda_true_std))
  P_true = Q_true%*%t(Q_true)
  
  n_ang = min(K_true, K_fit)
  angle_mat = matrix(NA, n_use, n_ang)
  proj_dist = numeric(n_use)
  lp = numeric(n_use)
  chain_id = integer(n_use)
  
  for (s in 1:n_use) {
    it = use_draws$iter[s]
    ch = use_draws$chain[s]
    
    L_s = matrix(as.numeric(lambda_arr[it, ch, Lambda_vars]), P, K_fit)
    psi_s = as.numeric(psi_arr[it, ch, psi_vars])
    sigb_s = as.numeric(sigb_arr[it, ch, sigb_vars])
    sige_s = as.numeric(sige_arr[it, ch, sige_vars])
    
    B_s = L_s%*%t(L_s)
    Sigma_a_s = B_s+diag(psi_s^2, P)
    
    Lambda_mean = Lambda_mean+L_s/n_use
    B_a_mean = B_a_mean+B_s/n_use
    Sigma_a_mean = Sigma_a_mean+Sigma_a_s/n_use
    R_a_mean = R_a_mean+cov2cor(Sigma_a_s)/n_use
    
    psi2_mean = psi2_mean+psi_s^2/n_use
    sigb2_mean = sigb2_mean+sigb_s^2/n_use
    sige2_mean = sige2_mean+sige_s^2/n_use
    pop_var_mean = pop_var_mean+(diag(Sigma_a_s)+sigb_s^2+sige_s^2)/n_use
    
    eig_common[s, ] = eigen(B_s, symmetric = TRUE, only.values = TRUE)$values[1:K_fit]
    
    Q_s = qr.Q(qr(L_s))
    sv = svd(t(Q_true)%*%Q_s)$d
    sv = pmin(pmax(sv, 0), 1)
    
    angle_mat[s, ] = acos(sv[1:n_ang])*180/pi
    
    P_s = Q_s%*%t(Q_s)
    proj_dist[s] = sqrt(sum((P_s-P_true)^2))
    lp[s] = as.numeric(lp_arr[it, ch, "lp__"])
    chain_id[s] = ch
  }
  
  off = upper.tri(R_a_true)
  
  max_rhat = max(diag_summ$rhat, na.rm = TRUE)
  min_ess = min(diag_summ$ess_bulk, na.rm = TRUE)
  
  lambda_rmse = NA
  lambda_cor = NA
  
  if (K_fit==K_true) {
    lambda_rmse = sqrt(mean((Lambda_mean-Lambda_true_std)^2))
    lambda_cor = cor(as.vector(Lambda_true_std), as.vector(Lambda_mean))
  }
  
  summary = data.frame(
    K_fit = K_fit,
    elapsed_sec = unname(time_fit["elapsed"]),
    max_rhat = max_rhat,
    min_ess = min_ess,
    cor_Ra = cor(R_a_true[off], R_a_mean[off]),
    rmse_Ra = sqrt(mean((R_a_mean[off]-R_a_true[off])^2)),
    cor_Ba = cor(B_a_true[off], B_a_mean[off]),
    rmse_Ba = sqrt(mean((B_a_mean[off]-B_a_true[off])^2)),
    cor_lambda = lambda_cor,
    rmse_lambda = lambda_rmse,
    rmse_common_diag = sqrt(mean((diag(B_a_mean)-diag(B_a_true))^2)),
    rmse_psi2 = sqrt(mean((psi2_mean-Psi_a_true)^2)),
    rmse_Sigma_a_diag = sqrt(mean((diag(Sigma_a_mean)-diag(Sigma_a_true))^2)),
    rmse_Db = sqrt(mean((sigb2_mean-Db_true)^2)),
    rmse_De = sqrt(mean((sige2_mean-De_true)^2)),
    rmse_pop_var = sqrt(mean((pop_var_mean-pop_var_true)^2)),
    median_proj_dist = median(proj_dist),
    median_angle1 = median(angle_mat[, 1]),
    median_angle2 = if (ncol(angle_mat)>=2) median(angle_mat[, 2]) else NA
  )
  
  list(
    K_fit = K_fit,
    fit = fit,
    summary = summary,
    R_a_mean = R_a_mean,
    B_a_mean = B_a_mean,
    Sigma_a_mean = Sigma_a_mean,
    eig_common = eig_common,
    angle_mat = angle_mat,
    proj_dist = proj_dist,
    chain_id = chain_id,
    lp = lp
  )
}

# -----------------------------
# Run all ranks
# -----------------------------

results = lapply(K_values, fit_one_rank)
names(results) = paste0("K", K_values)

summary_table = do.call(rbind, lapply(results, function(x) x$summary))
print(summary_table)

# -----------------------------
# Plots
# -----------------------------

par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

plot(summary_table$K_fit, summary_table$rmse_Ra,
     type = "b", pch = 19,
     xlab = "Fitted K",
     ylab = "RMSE",
     main = expression(cor(Sigma[a])))
abline(v = K_true, lty = 2)

plot(summary_table$K_fit, summary_table$rmse_Ba,
     type = "b", pch = 19,
     xlab = "Fitted K",
     ylab = "RMSE",
     main = expression(Lambda*Lambda^T))
abline(v = K_true, lty = 2)

plot(summary_table$K_fit, summary_table$elapsed_sec,
     type = "b", pch = 19,
     xlab = "Fitted K",
     ylab = "Seconds",
     main = "Runtime")
abline(v = K_true, lty = 2)

plot(summary_table$K_fit, summary_table$max_rhat,
     type = "b", pch = 19,
     xlab = "Fitted K",
     ylab = "Max Rhat",
     main = "Sampler diagnostic")
abline(h = 1.01, lty = 2)
abline(v = K_true, lty = 2)

# Correlation recovery by fitted K

par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))

off = upper.tri(R_a_true)

for (r in seq_along(results)) {
  plot(R_a_true[off], results[[r]]$R_a_mean[off],
       pch = 19,
       xlab = "True player-effect correlation",
       ylab = "Posterior mean",
       main = paste0("K = ", results[[r]]$K_fit))
  abline(0, 1, lty = 2)
}

# Common covariance eigenvalues

eig_true = eigen(B_a_true, symmetric = TRUE, only.values = TRUE)$values

ylim_eig = range(c(
  eig_true[1:8],
  unlist(lapply(results, function(x) {
    eigen(x$B_a_mean, symmetric = TRUE, only.values = TRUE)$values[1:8]
  }))
))

plot(eig_true[1:8],
     type = "b", pch = 19,
     ylim = ylim_eig,
     xlab = "Index",
     ylab = "Eigenvalue",
     main = expression("Eigenvalues of " * Lambda*Lambda^T))

for (r in seq_along(results)) {
  eig_r = eigen(results[[r]]$B_a_mean, symmetric = TRUE, only.values = TRUE)$values
  lines(eig_r[1:8], type = "b", pch = 19, lty = r+1)
}

legend("topright",
       legend = c("Truth", paste0("K=", K_values)),
       lty = c(1, 2:(length(K_values)+1)),
       pch = 19,
       bty = "n")

par(mfrow = c(1, 1))