data {
  int N;
  int T;
  int K;
  array[N*T] int ID;
  matrix[N*T, K] x;
  array[N*T] int seizures;
}

parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> u_sd;
  vector[N] u;
}

transformed parameters {
  vector[N*T] lambda = u[ID] + alpha + x * beta;
}

model {
  alpha ~ normal(0, 5);
  beta ~ std_normal();
  u_sd ~ cauchy(0, 5);
  u ~ normal(0, u_sd);
  seizures ~ poisson_log(lambda);
}

generated quantities {
  array[N*T] int ypred = poisson_log_rng(lambda);
}
