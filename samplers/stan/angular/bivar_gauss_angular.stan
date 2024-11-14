functions {
  #include gauge_fcns_euc.stanfunctions
  real est_vol(int n_vol, vector sum_term, vector sqrt_term, real pars) {
    vector[n_vol] gx = (sum_term - 2 * pars * sqrt_term) / (1 - pars^2);
    vector[n_vol] approx_indicator;
    real k = 15;
    for(i in 1:n_vol) {
     approx_indicator[i] = inv_logit(k * (1.0 - gx[i]));
    }
    return(mean(approx_indicator));
  }

  real angular_lpdf(array[] real angle, real pars, int dim, real L_volume) {
    int N = num_elements(angle);
    vector[N] angle_vec = to_vector(angle);
    return(-dim * sum(log(gauss_gauge(angle_vec, (1-angle_vec), pars))) - N * (log(dim) + log(L_volume)));
  }
}

data {
  int<lower=1> N;
  int<lower=1> n_grid;
  int<lower=1> d;
  array[N] real<lower=0, upper=1> w1;
  vector[n_grid] x1;
  vector[n_grid] x2;
}

transformed data {
  vector[n_grid] sum_x1_x2 = x1 + x2;
  vector[n_grid] sqrt_x1_x2 = sqrt(x1 .* x2);
}

parameters {
  real<lower=0, upper=1> dep;
}

transformed parameters {
  real<lower=0, upper=1> L = est_vol(n_grid, sum_x1_x2, sqrt_x1_x2, dep);
}

model {
  dep ~ uniform(0, 1);
  target += angular_lpdf(w1 | dep, d, L);
}

// generated quantities {
//   vector[N] log_lik;
//   for (n in 1:N) {
//     log_lik[n] = angular_lpdf(w1[n] | dep, d, L);
//   }
// }
