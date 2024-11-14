functions {
  #include gauge_fcns.stanfunctions
  #include likelihoodGamma.stanfunctions
  
  real angular_lpdf(real angle, real pars, int dim, real L_volume) {
        return(-dim * log(logistic_gauge(angle, pars)) - log(dim) - log(L_volume));
  }
}

data {
  int<lower=1> N;
  array[N] real<lower=0> R;
  array[N] real<lower=0, upper=1> W;
  array[N] real<lower=0> r0_w;
}

transformed data {
  int<lower=1> d = 2;
}

parameters {
  real<lower=0> alpha;
  real<lower=0, upper=1> dep1;
  real<lower=0, upper=1> dep2;
}

transformed parameters {
  // real<lower=0, upper=1> dep2 = dep1;
  real<lower=0, upper=1> L = dep1;
}

model {
  alpha ~ gamma(4, 2);
  dep1 ~ uniform(0, 1);
  dep2 ~ uniform(0, 1);
  
  for (n in 1:N) {
    target += (angular_lpdf(W[n] | dep1, d, L) + cens_gamma_lpdf(R[n] | r0_w[n], alpha, logistic_gauge(W[n], dep2)));
  }
}

// generated quantities {
//   vector[N] log_lik;
//   for (n in 1:N) {
//     log_lik[n] = angular_lpdf(w1[n] | dep, d, L);
//   }
// }
