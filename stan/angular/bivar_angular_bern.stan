data {
  int<lower=1> N;
  array[N] real<lower=0, upper=1> w1;
  int<lower=1> K;
}

parameters {
  simplex[K] weights; 
}

model {
  vector[K] log_weights = log(weights); 
  for (n in 1:N) {
    vector[K] lps = log_weights;
    for (j in 1:K) {
      real alpha = j;
      real beta = K - j + 1;
      lps[j] += beta_lpdf(w1[n] | alpha, beta);
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
