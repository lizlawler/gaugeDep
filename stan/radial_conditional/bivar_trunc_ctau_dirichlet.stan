functions {
  #include gauge_fcns.stanfunctions
  #include likelihoodGamma.stanfunctions
}

data {
  int<lower=1> N;
  array[N] real<lower=0> R;
  array[N] real<lower=0, upper=1> W;
  int<lower=1> n0_ctau;
  array[N] real<lower=0> r0_w_ctau;
  array[n0_ctau] int<lower=1> idx_ctau;
}

transformed data {
  array[n0_ctau] real<lower=0> R_trunc = R[idx_ctau];
  array[n0_ctau] real<lower=0, upper=1> W_trunc = W[idx_ctau];
  array[n0_ctau] real<lower=0> r0_w_trunc = r0_w_ctau[idx_ctau];
}

parameters {
  real<lower=0> alpha;
  real<lower=0> theta1;
  real<lower=0> theta2;
}

model {
  alpha ~ gamma(4, 2);
  theta1 ~ student_t(4, 0, 4);
  theta2 ~ student_t(4, 0, 2);
  
  for (n in 1:n0_ctau) {
    target += trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, dirichlet_gauge(W_trunc[n], theta1, theta2));
  }
}

generated quantities {
  vector[n0_ctau] log_lik;
  for (n in 1:n0_ctau) {
    log_lik[n] = trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, dirichlet_gauge(W_trunc[n], theta1, theta2));
  }
}
