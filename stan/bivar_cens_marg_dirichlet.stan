functions {
  #include gauge_fcns.stanfunctions
  #include likelihoodGamma.stanfunctions
}

data {
  int<lower=1> N;
  array[N] real<lower=0> R;
  array[N] real<lower=0, upper=1> W;
  array[N] real<lower=0> r0_w;
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
  
  for (n in 1:N) {
    target += cens_gamma_lpdf(R[n] | r0_w[n], alpha, dirichlet_gauge(W[n], theta1, theta2));
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = cens_gamma_lpdf(R[n] | r0_w[n], alpha, dirichlet_gauge(W[n], theta1, theta2));
  }
}
