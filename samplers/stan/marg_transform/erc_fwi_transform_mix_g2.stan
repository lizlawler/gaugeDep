functions {
  #include gpd_fcns.stanfunctions
  #include g2_fcns.stanfunctions
  #include trunc_expo_fcns.stanfunctions
  #include trunc_norm_fcns.stanfunctions
}

data {
  int<lower=1> N;
  int<lower=1> D;
  vector<lower=0>[N] erc;
  vector<lower=0>[N] fwi;
}

transformed data {
  // real<lower=0> fwi_trunc = 45;
  // real<lower=0> erc_trunc = 100;
  vector<lower=0>[D] max_vals;
  max_vals[1] = max(erc);
  max_vals[2] = max(fwi);
}

parameters {
  // egpd params (1 = ERC, 2 = FWI)
  vector[D] xi_prime;
  vector[D] kappa1_prime;
  vector[D] kappa2_prime;
  vector[D] sigma_prime;
  vector[D] prob_prime;
  
  // // exponential distribution params (for FWI)
  // real<lower=0> rate;
}

transformed parameters{
  vector<lower=-1>[D] xi = log1p_exp(xi_prime) - 1;
  vector<lower=0>[D] kappa1 = exp(kappa1_prime);
  vector<lower=0>[D] kappa2 = exp(kappa2_prime);
  vector<lower=0>[D] sigma;
  vector<lower=0, upper=1>[D] prob = inv_logit(prob_prime);
  // real<lower=0> sigma_norm = sqrt(sigma_sq);
  
  for(i in 1:D) {
    if(xi[i] < 0) {
      real min_sigma = max_vals[i] * (-xi[i]);
      sigma[i] = min_sigma + exp(sigma_prime[i]);
    } else {
      sigma[i] = exp(sigma_prime[i]);
    }
  }
  
  
}

model {
  to_vector(xi_prime) ~ std_normal();
  to_vector(kappa1_prime) ~ std_normal();
  to_vector(kappa2_prime) ~ std_normal();
  to_vector(sigma_prime) ~ normal(0, 4);
  
  // rate ~ lognormal(0,2);
  // pi_prob ~ normal(0, 4);
  to_vector(prob) ~ normal(0, 4);
  
  for (n in 1:N) {
    target += egpd_lpdf(erc[n] | sigma[1], xi[1], kappa1[1], kappa2[1], prob[1]);
    target += egpd_lpdf(fwi[n] | sigma[2], xi[2], kappa1[2], kappa2[2], prob[2]);
  }
}
