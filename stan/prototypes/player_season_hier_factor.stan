
data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> J;
  int<lower=1> T;
  int<lower=1> R_player;
  int<lower=1> R_season;
  array[N] int<lower=1, upper=J> player_id;
  array[N] int<lower=1, upper=T> season_id;
  matrix[N, D] Y;
}
parameters {
  vector[D] mu;

  matrix[J, R_player] a_raw;
  matrix[T, R_season] b_raw;

  matrix[D, R_player] Lambda_player;
  matrix[D, R_season] Lambda_season;

  vector<lower=0>[D] sigma_y;
}
transformed parameters {
  matrix[J, R_player] a;
  matrix[T, R_season] b;
  a = a_raw;
  b = b_raw;
}
model {
  // priors
  mu ~ normal(0, 1.5);
  to_vector(a_raw) ~ normal(0, 1);
  to_vector(b_raw) ~ normal(0, 1);
  to_vector(Lambda_player) ~ normal(0, 0.5);
  to_vector(Lambda_season) ~ normal(0, 0.5);
  sigma_y ~ normal(0, 0.5);

  // likelihood
  for (n in 1:N) {
    vector[D] mean_n;
    mean_n = mu + Lambda_player * to_vector(a[player_id[n]]') + Lambda_season * to_vector(b[season_id[n]]');
    Y[n] ~ normal(mean_n', sigma_y');
  }
}
generated quantities {
  matrix[N, D] Y_rep;
  vector[N] log_lik;

  for (n in 1:N) {
    vector[D] mean_n;
    mean_n = mu + Lambda_player * to_vector(a[player_id[n]]') + Lambda_season * to_vector(b[season_id[n]]');
    for (d in 1:D) {
      Y_rep[n, d] = normal_rng(mean_n[d], sigma_y[d]);
    }
    log_lik[n] = normal_lpdf(Y[n] | mean_n', sigma_y');
  }
}

