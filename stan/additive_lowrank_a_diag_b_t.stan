// Low-rank Sigma_a + diagonal Sigma_b with Student-t residuals.
//
// Identical to additive_lowrank_a_diag_b.stan except:
//   - No global mean mu (data is Z-scaled before fitting).
//   - Residual distribution is Student-t(nu) instead of Normal.
//     nu is estimated with prior Gamma(2, 0.1), as recommended by
//     Juarez & Steel (2010) and the Stan documentation.
//   - Scores are centred first (same numerical approach as the Normal model),
//     then A[I,P] and B[S,P] are precomputed as matrices so the model block
//     uses a single vectorised student_t call per observation instead of a
//     nested p-loop. This avoids inf−inf = NaN that arises when centering
//     large effects rather than small scores.
//
// If nu is large (> ~30) the t distribution is effectively Normal, so this
// model nests the normal-errors version. A small nu (< 10) indicates that
// heavier tails are needed (outlier player-seasons present).

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
  // Set to 1 to compute log_lik for LOO; 0 skips the N-loop in generated quantities
  // (saves ~N×P likelihood evaluations per draw on exploratory runs).
  int<lower=0, upper=1> compute_log_lik;
}

parameters {
  // Degrees of freedom for the t residuals
  real<lower=2> nu;

  // Player: low-rank factor scores + uniqueness residuals
  matrix[I, Q_a] eta_a_raw;
  matrix[I, P]   z_a_raw;

  // Player loading parameters
  vector<lower=0>[Q_a] lambda_a_diag;
  vector[N_lambda_a_free] lambda_a_free;
  vector<lower=0>[P] psi_a;

  // Season: plain diagonal
  matrix[S, P] z_b;
  vector<lower=0>[P] sigma_b;

  // Residual scale
  vector<lower=sigma_floor>[P] sigma_e;
}

transformed parameters {
  matrix[P, Q_a] Lambda_a = make_lower_tri_loadings(P, Q_a, lambda_a_diag, lambda_a_free);

  // Precompute player effect matrix A[I, P] and season effect matrix B[S, P].
  // Centering is applied to the RAW SCORES (not to the effects), which avoids
  // the inf - inf = NaN instability that occurs when centering large effect matrices.

  // 1. Centre factor scores: eta_c[i,q] = eta_a_raw[i,q] - mean_q(eta_a_raw)
  row_vector[Q_a] eta_a_bar = rep_row_vector(1.0 / I, I) * eta_a_raw;
  matrix[I, Q_a] eta_c = eta_a_raw - rep_matrix(eta_a_bar, I);

  // 2. Centre uniqueness scores: z_a_c[i,p] = z_a_raw[i,p] - mean_p(z_a_raw)
  row_vector[P] z_a_bar = rep_row_vector(1.0 / I, I) * z_a_raw;
  matrix[I, P] z_a_c = z_a_raw - rep_matrix(z_a_bar, I);

  // 3. Player effects: A = eta_c * Lambda_a' + z_a_c .* psi_a  (I × P matrix)
  matrix[I, P] A = eta_c * Lambda_a' + z_a_c .* rep_matrix(psi_a', I);

  // 4. Centre season scores: z_b_c[s,p] = z_b[s,p] - mean_p(z_b)
  row_vector[P] z_b_bar = rep_row_vector(1.0 / S, S) * z_b;
  matrix[S, P] z_b_c = z_b - rep_matrix(z_b_bar, S);

  // 5. Season effects: B = z_b_c .* sigma_b  (S × P matrix)
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

  // One vectorised student_t call per observation (over P features).
  for (n in 1:N) {
    vector[P] mean_n = to_vector(A[player_index[n]] + B[season_index[n]]);
    Y[n]' ~ student_t(nu, mean_n, sigma_e);
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
    // For t(nu), Var(ε) = sigma_e^2 * nu/(nu-2).
    var_e[p] = square(sigma_e[p]) * nu / (nu - 2);

    real total_var = Sigma_a_diag[p] + var_b[p] + var_e[p];
    prop_a[p] = Sigma_a_diag[p] / total_var;
    prop_b[p] = var_b[p]         / total_var;
    prop_e[p] = var_e[p]         / total_var;
  }

  if (compute_log_lik) {
    for (n in 1:N) {
      vector[P] mean_n = to_vector(A[player_index[n]] + B[season_index[n]]);
      log_lik[n] = student_t_lpdf(Y[n]' | nu, mean_n, sigma_e);
    }
  } else {
    log_lik = rep_vector(0, N);
  }
}
