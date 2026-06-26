data {
  int<lower=1> N_obs;
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> P;

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

  real<lower=0> scale_scale;
}

parameters {
  matrix[I, P] z_a;
  matrix[J, P] z_b;

  vector<lower=1e-6>[P] sigma_a;
  vector<lower=1e-6>[P] sigma_b;
  vector<lower=1e-6>[P] sigma_e;
}

model {
  to_vector(z_a) ~ std_normal();
  to_vector(z_b) ~ std_normal();

  sigma_a ~ student_t(3, 0, scale_scale);
  sigma_b ~ student_t(3, 0, scale_scale);
  sigma_e ~ student_t(3, 0, scale_scale);

  for (q in 1:N_train) {
    int n = train_obs[q];
    int p = train_var[q];
    real mu = sigma_a[p] * z_a[player[n], p]
      + sigma_b[p] * z_b[season[n], p];

    y_train[q] ~ normal(mu, sigma_e[p]);
  }
}

generated quantities {
  vector[N_test] log_lik_test;

  for (q in 1:N_test) {
    int n = test_obs[q];
    int p = test_var[q];
    real mu = sigma_a[p] * z_a[player[n], p]
      + sigma_b[p] * z_b[season[n], p];

    log_lik_test[q] = normal_lpdf(y_test[q] | mu, sigma_e[p]);
  }
}
