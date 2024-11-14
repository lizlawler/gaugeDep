data {
  int<lower=1> N;
  array[N] real<lower=0, upper=1> w1;
  // int<lower=1> K;
}

transformed data {
  int<lower=1> K = 3;
}

parameters {
  simplex[K] weights; 
  // ordered[K] mu;
  vector[K] mu_real;
  real<lower=0> beta_tau;
  vector<lower=0>[K] inv_tau;
}

transformed parameters {
  vector<lower=0>[K] tau = 1 / inv_tau;
  vector<lower=0,upper=1>[K] mu = inv_logit(mu_real);
  vector<lower=0>[K] alpha = mu .* tau;
  vector<lower=0>[K] beta = (1-mu) .* tau;
}

model {
  vector[K] log_weights = log(weights); 
  mu ~ uniform(0,1);
  beta_tau ~ exponential(1/8);
  inv_tau ~ gamma(2, beta_tau);
  
  for (n in 1:N) {
    vector[K] lps = log_weights;
    for (j in 1:K) {
      lps[j] += beta_lpdf(w1[n] | alpha[j], beta[j]);
    }
    target += log_sum_exp(lps);
  }
}

// generated quantities {
//   vector[N] log_lik;
//   for (n in 1:N) {
//     log_lik[n] = angular_lpdf(w1[n] | dep, d, L);
//   }
// }
