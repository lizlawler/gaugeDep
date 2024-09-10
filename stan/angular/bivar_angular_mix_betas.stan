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
  vector<lower=0, upper=1>[K] mu;
  vector<lower=0>[K] tau;
}

transformed parameters {
  vector<lower=0>[K] alpha = mu .* tau;
  vector<lower=0>[K] beta = (1-mu) .* tau;
}

model {
  vector[K] log_weights = log(weights); 
  mu ~ uniform(0,1);
  tau ~ lognormal(0, 2);
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
