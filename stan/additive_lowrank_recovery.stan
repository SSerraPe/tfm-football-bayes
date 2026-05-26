functions {
  matrix make_anchor_loadings(int P, int Q, array[] int anchor, vector positive_diag, vector free_values) {
    matrix[P, Q] L = rep_matrix(0, P, Q);
    int pos = 1;

    for (q in 1:Q) {
      for (p in 1:P) {
        int is_anchor_row = 0;
        for (k in 1:Q) {
          if (p == anchor[k]) {
            is_anchor_row = 1;
          }
        }

        if (p == anchor[q]) {
          L[p, q] = positive_diag[q];
        } else if (is_anchor_row == 0) {
          L[p, q] = free_values[pos];
          pos += 1;
        }
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
  int<lower=1> Q_b;
  int<lower=0> N_lambda_a_free;
  int<lower=0> N_lambda_b_free;
  matrix[N, P] Y;
  array[N] int<lower=1, upper=I> player_index;
  array[N] int<lower=1, upper=S> season_index;
  array[Q_a] int<lower=1, upper=P> lambda_a_anchor;
  array[Q_b] int<lower=1, upper=P> lambda_b_anchor;
  real<lower=0> sigma_floor;
}

parameters {
  vector[P] mu;
  matrix[I, Q_a] eta_a_raw;
  matrix[S, Q_b] eta_b_raw;
  matrix[I, P] z_a_raw;
  matrix[S, P] z_b_raw;
  vector<lower=0>[Q_a] lambda_a_diag;
  vector[N_lambda_a_free] lambda_a_free;
  vector<lower=0>[Q_b] lambda_b_diag;
  vector[N_lambda_b_free] lambda_b_free;
  vector<lower=0>[P] psi_a;
  vector<lower=0>[P] psi_b;
  vector<lower=sigma_floor>[P] sigma_e;
}

transformed parameters {
  matrix[P, Q_a] Lambda_a = make_anchor_loadings(P, Q_a, lambda_a_anchor, lambda_a_diag, lambda_a_free);
  matrix[P, Q_b] Lambda_b = make_anchor_loadings(P, Q_b, lambda_b_anchor, lambda_b_diag, lambda_b_free);
  vector[Q_a] eta_a_bar;
  vector[Q_b] eta_b_bar;
  vector[P] z_a_bar;
  vector[P] z_b_bar;

  for (q in 1:Q_a) {
    eta_a_bar[q] = 0;
    for (i in 1:I) {
      eta_a_bar[q] += eta_a_raw[i, q];
    }
    eta_a_bar[q] /= I;
  }

  for (q in 1:Q_b) {
    eta_b_bar[q] = 0;
    for (s in 1:S) {
      eta_b_bar[q] += eta_b_raw[s, q];
    }
    eta_b_bar[q] /= S;
  }

  for (p in 1:P) {
    z_a_bar[p] = 0;
    z_b_bar[p] = 0;
    for (i in 1:I) {
      z_a_bar[p] += z_a_raw[i, p];
    }
    z_a_bar[p] /= I;
    for (s in 1:S) {
      z_b_bar[p] += z_b_raw[s, p];
    }
    z_b_bar[p] /= S;
  }
}

model {
  mu ~ normal(0, 1);
  to_vector(eta_a_raw) ~ std_normal();
  to_vector(eta_b_raw) ~ std_normal();
  to_vector(z_a_raw) ~ std_normal();
  to_vector(z_b_raw) ~ std_normal();
  lambda_a_diag ~ lognormal(log(0.6), 0.35);
  lambda_a_free ~ normal(0, 0.35);
  lambda_b_diag ~ lognormal(log(0.4), 0.35);
  lambda_b_free ~ normal(0, 0.25);
  psi_a ~ normal(0, 0.4);
  psi_b ~ normal(0, 0.3);
  sigma_e ~ normal(0, 0.7);

  for (n in 1:N) {
    for (p in 1:P) {
      real mean_np = mu[p]
        + psi_a[p] * (z_a_raw[player_index[n], p] - z_a_bar[p])
        + psi_b[p] * (z_b_raw[season_index[n], p] - z_b_bar[p]);
      for (q in 1:Q_a) {
        mean_np += Lambda_a[p, q] * (eta_a_raw[player_index[n], q] - eta_a_bar[q]);
      }
      for (q in 1:Q_b) {
        mean_np += Lambda_b[p, q] * (eta_b_raw[season_index[n], q] - eta_b_bar[q]);
      }
      Y[n, p] ~ normal(
        mean_np,
        sigma_e[p]
      );
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[P] Sigma_a_diag;
  vector[P] Sigma_b_diag;
  vector[P] Sigma_epsilon_diag;

  for (p in 1:P) {
    real factor_var_a = 0;
    real factor_var_b = 0;
    for (q in 1:Q_a) {
      factor_var_a += square(Lambda_a[p, q]);
    }
    for (q in 1:Q_b) {
      factor_var_b += square(Lambda_b[p, q]);
    }
    Sigma_a_diag[p] = factor_var_a + square(psi_a[p]);
    Sigma_b_diag[p] = factor_var_b + square(psi_b[p]);
    Sigma_epsilon_diag[p] = square(sigma_e[p]);
  }

  for (n in 1:N) {
    real lp = 0;
    for (p in 1:P) {
      real mean_np = mu[p]
        + psi_a[p] * (z_a_raw[player_index[n], p] - z_a_bar[p])
        + psi_b[p] * (z_b_raw[season_index[n], p] - z_b_bar[p]);
      for (q in 1:Q_a) {
        mean_np += Lambda_a[p, q] * (eta_a_raw[player_index[n], q] - eta_a_bar[q]);
      }
      for (q in 1:Q_b) {
        mean_np += Lambda_b[p, q] * (eta_b_raw[season_index[n], q] - eta_b_bar[q]);
      }
      lp += normal_lpdf(
        Y[n, p] |
        mean_np,
        sigma_e[p]
      );
    }
    log_lik[n] = lp;
  }
}
