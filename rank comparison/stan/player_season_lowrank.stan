data {
  int<lower=1> N_obs;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> P;
  int<lower=1, upper=P> K;

  array[N_obs] int<lower=1, upper=I> player;
  array[N_obs] int<lower=1, upper=J> season;

  int<lower=1> N_train;
  array[N_train] int<lower=1, upper=N_obs> train_obs;
  array[N_train] int<lower=1, upper=P> train_var;
  vector[N_train] y_train;

  int<lower=0> N_test;
  array[N_test] int<lower=1, upper=N_obs> test_obs;
  array[N_test] int<lower=1, upper=P> test_var;
  vector[N_test] y_test;

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
  matrix[P, K] Lambda = rep_matrix(0, P, K);
  vector<lower=0>[K] lambda_col_norm;

  {
    int pos = 1;

    for (h in 1:K) {
      Lambda[h, h] = lambda_diag[h];

      if (h < P) {
        for (p in (h + 1):P) {
          Lambda[p, h] = lambda_lower[pos];
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
  matrix[I, P] A_low = F * Lambda';

  to_vector(F) ~ std_normal();
  to_vector(z_u) ~ std_normal();
  to_vector(z_b) ~ std_normal();

  lambda_diag ~ normal(0, lambda_scale);
  lambda_lower ~ normal(0, lambda_scale);

  psi_a ~ student_t(3, 0, scale_scale);
  sigma_b ~ student_t(3, 0, scale_scale);
  sigma_e ~ student_t(3, 0, scale_scale);

  for (q in 1:N_train) {
    int n = train_obs[q];
    int p = train_var[q];
    real mu = A_low[player[n], p]
      + psi_a[p] * z_u[player[n], p]
      + sigma_b[p] * z_b[season[n], p];

    y_train[q] ~ normal(mu, sigma_e[p]);
  }
}

generated quantities {
  matrix[I, P] A_low = F * Lambda';
  vector[N_test] log_lik_test;

  for (q in 1:N_test) {
    int n = test_obs[q];
    int p = test_var[q];
    real mu = A_low[player[n], p]
      + psi_a[p] * z_u[player[n], p]
      + sigma_b[p] * z_b[season[n], p];

    log_lik_test[q] = normal_lpdf(y_test[q] | mu, sigma_e[p]);
  }
}
