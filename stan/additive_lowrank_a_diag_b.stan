functions {
  // Lower-triangular identification: L[q,q] > 0, L[p,q] free for p > q, L[p,q] = 0 for p < q.
  // Free values are stored column-major: col 1 rows 2..P, then col 2 rows 3..P, etc.
  matrix make_lower_tri_loadings(int P, int Q, vector positive_diag, vector free_values) {
    matrix[P, Q] L = rep_matrix(0, P, Q);
    int pos = 1;
    for (q in 1:Q) {
      L[q, q] = positive_diag[q];
      for (p in (q + 1):P) {
        L[p, q] = free_values[pos];
        pos += 1;
      }
    }
    return L;
  }
}

data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> I;
  int<lower=1> S;
  int<lower=1> Q_a;
  int<lower=0> N_lambda_a_free;
  matrix[N, P] Y;
  array[N] int<lower=1, upper=I> player_index;
  array[N] int<lower=1, upper=S> season_index;
  real<lower=0> sigma_floor;
}

parameters {
  vector[P] mu;

  // Player: low-rank factor scores + uniqueness residuals
  matrix[I, Q_a] eta_a_raw;
  matrix[I, P] z_a_raw;

  // Player loading parameters
  vector<lower=0>[Q_a] lambda_a_diag;
  vector[N_lambda_a_free] lambda_a_free;
  vector<lower=0>[P] psi_a;

  // Season: plain diagonal (as in additive_diagonal.stan)
  matrix[S, P] z_b;
  vector<lower=0>[P] sigma_b;

  // Residual
  vector<lower=sigma_floor>[P] sigma_e;
}

transformed parameters {
  matrix[P, Q_a] Lambda_a = make_lower_tri_loadings(P, Q_a, lambda_a_diag, lambda_a_free);

  // Column-centre all random effects for identifiability
  vector[Q_a] eta_a_bar;
  vector[P]   z_a_bar;
  vector[P]   z_b_bar;

  for (q in 1:Q_a) {
    eta_a_bar[q] = 0;
    for (i in 1:I) eta_a_bar[q] += eta_a_raw[i, q];
    eta_a_bar[q] /= I;
  }

  for (p in 1:P) {
    z_a_bar[p] = 0;
    z_b_bar[p] = 0;
    for (i in 1:I) z_a_bar[p] += z_a_raw[i, p];
    z_a_bar[p] /= I;
    for (s in 1:S) z_b_bar[p] += z_b[s, p];
    z_b_bar[p] /= S;
  }
}

model {
  mu             ~ normal(0, 1);
  to_vector(eta_a_raw) ~ std_normal();
  to_vector(z_a_raw)   ~ std_normal();
  to_vector(z_b)       ~ std_normal();

  lambda_a_diag ~ lognormal(log(0.7), 0.20);
  lambda_a_free ~ normal(0, 0.35);
  psi_a         ~ normal(0, 0.4);
  sigma_b       ~ normal(0, 0.5);
  sigma_e       ~ normal(0, 0.7);

  for (n in 1:N) {
    for (p in 1:P) {
      real mean_np = mu[p]
        + psi_a[p] * (z_a_raw[player_index[n], p] - z_a_bar[p])
        + sigma_b[p] * (z_b[season_index[n], p] - z_b_bar[p]);
      for (q in 1:Q_a) {
        mean_np += Lambda_a[p, q] * (eta_a_raw[player_index[n], q] - eta_a_bar[q]);
      }
      Y[n, p] ~ normal(mean_np, sigma_e[p]);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[P] Sigma_a_diag;
  vector[P] var_b;
  vector[P] var_e;
  vector[P] prop_a;
  vector[P] prop_b;
  vector[P] prop_e;

  for (p in 1:P) {
    real factor_var_a = 0;
    for (q in 1:Q_a) factor_var_a += square(Lambda_a[p, q]);
    Sigma_a_diag[p] = factor_var_a + square(psi_a[p]);
    var_b[p] = square(sigma_b[p]);
    var_e[p] = square(sigma_e[p]);

    real total_var = Sigma_a_diag[p] + var_b[p] + var_e[p];
    prop_a[p] = Sigma_a_diag[p] / total_var;
    prop_b[p] = var_b[p]         / total_var;
    prop_e[p] = var_e[p]         / total_var;
  }

  for (n in 1:N) {
    real lp = 0;
    for (p in 1:P) {
      real mean_np = mu[p]
        + psi_a[p] * (z_a_raw[player_index[n], p] - z_a_bar[p])
        + sigma_b[p] * (z_b[season_index[n], p] - z_b_bar[p]);
      for (q in 1:Q_a) {
        mean_np += Lambda_a[p, q] * (eta_a_raw[player_index[n], q] - eta_a_bar[q]);
      }
      lp += normal_lpdf(Y[n, p] | mean_np, sigma_e[p]);
    }
    log_lik[n] = lp;
  }
}
