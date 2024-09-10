// shared functions, data input values and shared parameter declaration
functions {
  #include gauge_fcns.stanfunctions
  #include likelihoodGamma.stanfunctions
  
  real est_vol(int n_vol, vector sum_term, vector sqrt_term, real pars) {
    vector[n_vol] gx = (sum_term - 2 * pars * sqrt_term) / (1 - pars^2);
    vector[n_vol] approx_indicator;
    real k = 15;
    for(i in 1:n_vol) {
     approx_indicator[i] = inv_logit(k * (1.0 - gx[i]));
    }
    return(mean(approx_indicator));
  }

  real angular_lpdf(real angle, real pars, int dim, real L_volume) {
    return(-dim * log(gauss_gauge(angle, pars)) - (log(dim) + log(L_volume)));
  }
}

data {
  int<lower=1> N;
  array[N] real<lower=0> R;
  array[N] real<lower=0, upper=1> W;
  array[N] real<lower=0> r0_w;
  int<lower=1> n_grid;
  vector[n_grid] x1;
  vector[n_grid] x2;
}

transformed data {
  int<lower=1> d = 2;
  vector[n_grid] sum_x1_x2 = x1 + x2;
  vector[n_grid] sqrt_x1_x2 = sqrt(x1 .* x2);
}

parameters {
  real<lower=0> alpha;
  real<lower=0, upper =1> dep1;
  real<lower=0, upper =1> dep2;
}

transformed parameters {
  // real<lower=0, upper =1> dep2 = dep1;
  real<lower=0, upper=1> L = est_vol(n_grid, sum_x1_x2, sqrt_x1_x2, dep1);
}

model {
  alpha ~ gamma(4, 2);
  dep1 ~ uniform(0, 1);
  dep2 ~ uniform(0, 1);
  
  for (n in 1:N) {
    target += (angular_lpdf(W[n] | dep1, d, L) + cens_gamma_lpdf(R[n] | r0_w[n], alpha, gauss_gauge(W[n], dep2)));
  }
}

// generated quantities {
//   vector[N] log_lik;
//   for (n in 1:N) {
//     log_lik[n] = cens_gamma_lpdf(R[n] | r0_w[n], alpha, gauss_gauge(W[n], dep));
//   }
// }
