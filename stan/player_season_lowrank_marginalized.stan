// Marginalized player-effects model (Normal residuals).
//
// Additive decomposition: y_n = a_{i[n]} + b_{j[n]} + e_n
//   a_i ~ N(0, Sigma_a = Lambda_a Lambda_a' + diag(psi_a^2))   [marginalized out]
//   b_j ~ N(0, diag(sigma_b^2))                                  [explicit via z_b]
//   e_n ~ N(0, diag(sigma_e^2))
//
// Marginalizing a_i, conditional on B, player i's residuals y_tilde_{i,s} = y_{i,s} - b_{j[n]}
// follow:
//   Cov(y_tilde_i) = J_{S_i} ox Sigma_a + I_{S_i} ox D_e     (D_e = diag(sigma_e^2))
//
// Woodbury identity gives the inverse and log-determinant in terms of the P x P matrix:
//   M_i = Sigma_a^{-1} + S_i * D_e^{-1}
//
// This reduces parameter dimensionality from ~I*P + I*K + S*P to ~K*P + 3*P + S*P,
// a ~100x reduction for typical values (I=1529, K=2, P=48, S=12).

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
  array[I] int<lower=1> S_count;     // number of player-season rows per player
  real<lower=0> sigma_floor;
}

parameters {
  vector<lower=0>[Q_a]        lambda_a_diag;
  vector[N_lambda_a_free]     lambda_a_free;
  vector<lower=0>[P]          psi_a;

  matrix[S, P]                z_b;
  vector<lower=0>[P]          sigma_b;

  vector<lower=sigma_floor>[P] sigma_e;
}

transformed parameters {
  matrix[P, Q_a] Lambda_a = make_lower_tri_loadings(P, Q_a, lambda_a_diag, lambda_a_free);

  // Centred season effects (same convention as additive_lowrank_a_diag_b.stan)
  row_vector[P] z_b_bar = rep_row_vector(1.0 / S, S) * z_b;
  matrix[S, P]  z_b_c   = z_b - rep_matrix(z_b_bar, S);
  matrix[S, P]  B        = z_b_c .* rep_matrix(sigma_b', S);
}

model {
  lambda_a_diag  ~ lognormal(log(0.7), 0.20);
  lambda_a_free  ~ normal(0, 0.35);
  psi_a          ~ normal(0, 0.4);
  to_vector(z_b) ~ std_normal();
  sigma_b        ~ normal(0, 0.5);
  sigma_e        ~ normal(0, 0.7);

  // Build Sigma_a and precompute stable quantities (O(P^3) once per iteration)
  matrix[P, P] Sigma_a     = tcrossprod(Lambda_a) + diag_matrix(square(psi_a));
  matrix[P, P] L_Sig       = cholesky_decompose(Sigma_a);
  real         log_det_Sig = 2.0 * sum(log(diagonal(L_Sig)));

  // Sigma_a^{-1}: invert via Cholesky factors
  matrix[P, P] L_inv       = mdivide_left_tri_low(L_Sig, diag_matrix(rep_vector(1.0, P)));
  matrix[P, P] Sigma_a_inv = L_inv' * L_inv;

  // Residual precision
  vector[P] sigma_e2     = square(sigma_e);
  vector[P] sigma_e2_inv = 1.0 ./ sigma_e2;
  real      log_det_De   = sum(log(sigma_e2));

  // Accumulate per-player centred residual statistics: B-subtracted sums and sumsq
  // y_tilde_sum[i, p] = sum_s (y_{i,s,p} - B[j[n], p])
  // y_tilde_ss[i, p]  = sum_s (y_{i,s,p} - B[j[n], p])^2
  matrix[I, P] y_tilde_sum = rep_matrix(0.0, I, P);
  matrix[I, P] y_tilde_ss  = rep_matrix(0.0, I, P);
  for (n in 1:N) {
    row_vector[P] r        = Y[n] - B[season_index[n]];
    y_tilde_sum[player_index[n]] += r;
    y_tilde_ss[player_index[n]]  += square(r);
  }

  // Per-player marginal log-likelihood contribution via Woodbury
  for (i in 1:I) {
    int Si = S_count[i];

    // M_i = Sigma_a^{-1} + S_i * diag(sigma_e^{-2})    [P x P, SPD]
    matrix[P, P] M_i = Sigma_a_inv + Si * diag_matrix(sigma_e2_inv);

    // log |Cov_i| = S_i * log|D_e| + log|Sigma_a| + log|M_i|
    real log_det_cov_i = Si * log_det_De + log_det_Sig + log_determinant(M_i);

    // Quadratic form:
    //   q1 = sum_p sigma_e^{-2}_p * y_tilde_ss[i,p]    (sum_s ||y_tilde_{i,s}||^2_{De^{-1}})
    //   q2 = y_tilde_sum_i' De^{-1} M_i^{-1} De^{-1} y_tilde_sum_i
    vector[P] yts_i      = (y_tilde_sum[i])';
    vector[P] De_yts     = sigma_e2_inv .* yts_i;
    real      q1         = dot_product(sigma_e2_inv, (y_tilde_ss[i])');
    vector[P] sol        = mdivide_left_spd(M_i, De_yts);
    real      q2         = dot_product(De_yts, sol);

    target += -0.5 * (Si * P * log(2.0 * pi()) + log_det_cov_i + q1 - q2);
  }
}

generated quantities {
  vector[P] Sigma_a_diag;
  vector[P] var_b;
  vector[P] var_e;
  vector[P] prop_a;
  vector[P] prop_b;
  vector[P] prop_e;

  for (p in 1:P) {
    Sigma_a_diag[p] = sum(square(Lambda_a[p])) + square(psi_a[p]);
    var_b[p]  = square(sigma_b[p]);
    var_e[p]  = square(sigma_e[p]);
    real tv   = Sigma_a_diag[p] + var_b[p] + var_e[p];
    prop_a[p] = Sigma_a_diag[p] / tv;
    prop_b[p] = var_b[p] / tv;
    prop_e[p] = var_e[p] / tv;
  }
}
