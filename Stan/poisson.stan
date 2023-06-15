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
}

transformed parameters {
  vector[N*T] lambda = alpha + x * beta;
}

model {
  alpha ~ normal(0, 5);
  beta ~ std_normal();
  seizures ~ poisson_log(lambda);
}

generated quantities {
  array[N*T] int ypred = poisson_log_rng(lambda);
}
