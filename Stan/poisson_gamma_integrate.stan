functions {
  real pois_gamma_integrand(real u, real unused, array[] real params,
                              array[] real unused2, array[] int y_i) {
    real lambda = params[1];
    real u_shape = params[2];
    real marg_lpdf = poisson_log_lpmf(y_i[1] | lambda + log(u))
                      + gamma_lpdf(u | u_shape, u_shape);
    return exp(marg_lpdf);
  }

  real poisson_gamma_lpmf(int y, real lambda, real u_shape) {
    real pdf = integrate_1d(pois_gamma_integrand,
                              0, positive_infinity(),
                              {lambda, u_shape},
                              {0},
                              {y}
                              );
    return log(pdf);
  }
}

data {
  int N;
  int T;
  array[N*T] int ID;
  vector[N*T] time;
  vector[N*T] baseline;
  vector[N*T] age;
  vector[N*T] treatment;
  array[N*T] int seizures;
}

transformed data {
 vector[N*T] base_x_treat = baseline .* treatment;
}

parameters {
  real alpha;
  real beta_treat;
  real beta_age;
  real beta_baseline;
  real beta_base_x_treat;
  real<lower=0> u_shape;
  vector<lower=0>[N] u;
}

transformed parameters {
  vector[N*T] lambda;

  lambda = alpha
            + beta_treat * treatment
            + beta_age * age
            + beta_baseline * baseline
            + beta_base_x_treat * base_x_treat;
}

model {
  alpha ~ normal(0, 5);
  beta_treat ~ std_normal();
  beta_age ~ std_normal();
  beta_baseline ~ std_normal();
  beta_base_x_treat ~ std_normal();
  u_shape ~ cauchy(0, 5);
  u ~ gamma(u_shape, u_shape);
  seizures ~ poisson_log(log(u[ID]) + lambda);
}

generated quantities {
  array[N*T] int ypred = poisson_log_rng(log(u[ID]) + lambda);
  vector[N*T] log_lik;

  for (i in 1:(N*T)) {
    log_lik[i] = poisson_gamma_lpmf(seizures[i] | lambda[i], u_shape);
  }
}
