// Low-rank Sigma_a + diagonal Sigma_b, Student-t residuals, minutes-scaled sigma_e.
//
// Extension of additive_lowrank_a_diag_b_t.stan with observation-level residual scale:
//   epsilon_{n,p} ~ t(nu, 0, sigma_e[p] * sqrt(m_ref / minutes_n))
//
// Motivation: per-90 features are ratio estimators (count / minutes * 90) whose
// sampling precision grows with playing time. For a Poisson event count X_n with
// rate lambda, Var(X_n / minutes_n) = lambda / minutes_n, so SD ~ 1/sqrt(minutes_n).
// The factor sqrt(m_ref / minutes_n) re-scales sigma_e to a reference playing time
// m_ref (typically the sample median), making sigma_e[p] the residual SD at m_ref.
//
// Fixed phi = 0.5 variant (see additive_lowrank_a_diag_b_t_mv_phi.stan for estimated phi).
// ICC is computed at m_ref: var_e[p] = sigma_e[p]^2 * nu/(nu-2), same formula as the
// base model — sigma_e[p] now represents the reference-scale residual.

functions {
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
  int<lower=0, upper=1> compute_log_lik;
  vector<lower=0>[N] minutes;   // playing time for each player-season
  real<lower=0> m_ref;          // reference minutes (e.g. sample median)
}

transformed data {
  // scale_n = sqrt(m_ref / minutes_n): >1 for low-minute, <1 for high-minute
  vector[N] sigma_scale = sqrt(m_ref ./ minutes);
}

parameters {
  real<lower=2> nu;

  matrix[I, Q_a] eta_a_raw;
  matrix[I, P]   z_a_raw;

  vector<lower=0>[Q_a] lambda_a_diag;
  vector[N_lambda_a_free] lambda_a_free;
  vector<lower=0>[P] psi_a;

  matrix[S, P] z_b;
  vector<lower=0>[P] sigma_b;

  vector<lower=sigma_floor>[P] sigma_e;
}

transformed parameters {
  matrix[P, Q_a] Lambda_a = make_lower_tri_loadings(P, Q_a, lambda_a_diag, lambda_a_free);

  row_vector[Q_a] eta_a_bar = rep_row_vector(1.0 / I, I) * eta_a_raw;
  matrix[I, Q_a] eta_c = eta_a_raw - rep_matrix(eta_a_bar, I);

  row_vector[P] z_a_bar = rep_row_vector(1.0 / I, I) * z_a_raw;
  matrix[I, P] z_a_c = z_a_raw - rep_matrix(z_a_bar, I);

  matrix[I, P] A = eta_c * Lambda_a' + z_a_c .* rep_matrix(psi_a', I);

  row_vector[P] z_b_bar = rep_row_vector(1.0 / S, S) * z_b;
  matrix[S, P] z_b_c = z_b - rep_matrix(z_b_bar, S);

  matrix[S, P] B = z_b_c .* rep_matrix(sigma_b', S);
}

model {
  nu ~ gamma(2, 0.1);

  to_vector(eta_a_raw) ~ std_normal();
  to_vector(z_a_raw)   ~ std_normal();
  to_vector(z_b)       ~ std_normal();

  lambda_a_diag ~ lognormal(log(0.7), 0.20);
  lambda_a_free ~ normal(0, 0.35);
  psi_a         ~ normal(0, 0.4);
  sigma_b       ~ normal(0, 0.5);
  sigma_e       ~ normal(0, 0.7);

  for (n in 1:N) {
    vector[P] mean_n  = to_vector(A[player_index[n]] + B[season_index[n]]);
    vector[P] scale_n = sigma_e * sigma_scale[n];
    Y[n]' ~ student_t(nu, mean_n, scale_n);
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
    Sigma_a_diag[p] = sum(square(Lambda_a[p])) + square(psi_a[p]);
    var_b[p] = square(sigma_b[p]);
    // var_e at the reference scale m_ref (phi=0.5, scale=1 at m_ref)
    var_e[p] = square(sigma_e[p]) * nu / (nu - 2);

    real total_var = Sigma_a_diag[p] + var_b[p] + var_e[p];
    prop_a[p] = Sigma_a_diag[p] / total_var;
    prop_b[p] = var_b[p]         / total_var;
    prop_e[p] = var_e[p]         / total_var;
  }

  if (compute_log_lik) {
    for (n in 1:N) {
      vector[P] mean_n  = to_vector(A[player_index[n]] + B[season_index[n]]);
      vector[P] scale_n = sigma_e * sigma_scale[n];
      log_lik[n] = student_t_lpdf(Y[n]' | nu, mean_n, scale_n);
    }
  } else {
    log_lik = rep_vector(0, N);
  }
}
