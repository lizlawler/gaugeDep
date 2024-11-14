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

transformed data {
  int<lower=1> K = 3;
}

parameters {
  // Gamma censored likelihood parameters
  real<lower=0> alpha;
  real<lower=0, upper =1> dep;
  
  // Mixture of beta density location and scale parameters
  simplex[K] weights; 
  vector<lower=0, upper=1>[K] mu;
  vector<lower=0>[K] tau;
}

transformed parameters {
  // Mixture of beta density true parameters
  vector<lower=0>[K] alpha_angles = mu .* tau;
  vector<lower=0>[K] beta_angles = (1-mu) .* tau;
}

model {
  alpha ~ gamma(4, 2);
  dep ~ uniform(0, 1);
  mu ~ uniform(0,1);
  tau ~ lognormal(0, 2);
  
  vector[K] log_weights = log(weights); 
  for (n in 1:N) {
    vector[K] lps = log_weights;
    for (j in 1:K) {
      lps[j] += beta_lpdf(W[n] | alpha_angles[j], beta_angles[j]);
    }
    target += log_sum_exp(lps);
    target += cens_gamma_lpdf(R[n] | r0_w[n], alpha, logistic_gauge(W[n], dep));
  }
}
