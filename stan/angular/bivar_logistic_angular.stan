functions {
  #include gauge_fcns_euc.stanfunctions
  // real est_vol_rng(int n1_vol, real pars) {
  //   vector[n1_vol] alpha = rep_vector(0, n1_vol);
  //   vector[n1_vol] beta = rep_vector(1, n1_vol);
  //   array[n1_vol] real x1 = uniform_rng(alpha, beta);
  //   array[n1_vol] real x2 = uniform_rng(alpha, beta);
  //   array[n1_vol] real gx;
  //   for(i in 1:n1_vol) {
  //    gx[i] = (logistic_gauge(x1[i], x2[i], pars) <= 1.0) ? 1 : 0;
  //   }
  //   return(sum(gx)/n1_vol);
  // }
  // 
  // real mc_volume(int n_mc_vol, int n1_vol, real pars) {
  //   vector[n_mc_vol] est_volume_n1;
  //   for(n in 1:n_mc_vol) {
  //     est_volume_n1[n] = est_vol_rng(n1_vol, pars);
  //   }
  //   return(mean(est_volume_n1));
  // }
  
  real angular_lpdf(real angle, real pars, int dim, real L_volume) {
        return(-dim * log(logistic_gauge(angle, (1-angle), pars)) - log(dim) - log(L_volume));
  }
}

data {
  int<lower=1> N;
  int<lower=1> d;
  array[N] real<lower=0, upper=1> w1;
  // real<lower=0, upper=1> L;
}

// transformed data {
//   int n1_vol = 500;
//   int n_mc_vol = 10000;
// }

parameters {
  real<lower=0, upper=1> dep;
}

transformed parameters {
  real<lower=0, upper=1> L = dep;
}

model {
  dep ~ uniform(0, 1);
  
  for (n in 1:N) {
    target += angular_lpdf(w1[n] | dep, d, L);
  }
}

// generated quantities {
//   vector[N] log_lik;
//   for (n in 1:N) {
//     log_lik[n] = angular_lpdf(w1[n] | dep, d, L);
//   }
// }
