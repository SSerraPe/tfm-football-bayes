data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> I;
  int<lower=1> S;
  matrix[N, P] Y;
  array[N] int<lower=1, upper=I> player_index;
  array[N] int<lower=1, upper=S> season_index;
  real<lower=0> sigma_floor;
}

parameters {
  vector[P] mu;
  matrix[I, P] z_a;
  matrix[S, P] z_b;
  vector<lower=0>[P] sigma_a;
  vector<lower=0>[P] sigma_b;
  vector<lower=sigma_floor>[P] sigma_e;
}

transformed parameters {
  vector[P] z_a_bar;
  vector[P] z_b_bar;

  for (p in 1:P) {
    z_a_bar[p] = 0;
    z_b_bar[p] = 0;
    for (i in 1:I) {
      z_a_bar[p] += z_a[i, p];
    }
    z_a_bar[p] /= I;
    for (s in 1:S) {
      z_b_bar[p] += z_b[s, p];
    }
    z_b_bar[p] /= S;
  }
}

model {
  mu ~ normal(0, 1);
  to_vector(z_a) ~ std_normal();
  to_vector(z_b) ~ std_normal();
  sigma_a ~ normal(0, 0.7);
  sigma_b ~ normal(0, 0.5);
  sigma_e ~ normal(0, 0.7);

  for (n in 1:N) {
    for (p in 1:P) {
      Y[n, p] ~ normal(
        mu[p]
          + sigma_a[p] * (z_a[player_index[n], p] - z_a_bar[p])
          + sigma_b[p] * (z_b[season_index[n], p] - z_b_bar[p]),
        sigma_e[p]
      );
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[P] var_a;
  vector[P] var_b;
  vector[P] var_e;
  vector[P] prop_a;
  vector[P] prop_b;
  vector[P] prop_e;

  for (p in 1:P) {
    real total_var;
    var_a[p] = square(sigma_a[p]);
    var_b[p] = square(sigma_b[p]);
    var_e[p] = square(sigma_e[p]);
    total_var = var_a[p] + var_b[p] + var_e[p];
    prop_a[p] = var_a[p] / total_var;
    prop_b[p] = var_b[p] / total_var;
    prop_e[p] = var_e[p] / total_var;
  }

  for (n in 1:N) {
    real lp = 0;
    for (p in 1:P) {
      lp += normal_lpdf(
        Y[n, p] |
        mu[p]
          + sigma_a[p] * (z_a[player_index[n], p] - z_a_bar[p])
          + sigma_b[p] * (z_b[season_index[n], p] - z_b_bar[p]),
        sigma_e[p]
      );
    }
    log_lik[n] = lp;
  }
}
