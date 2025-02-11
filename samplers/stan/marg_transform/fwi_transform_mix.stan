functions {
  #include gpd_fcns.stanfunctions
  #include g1_fcns.stanfunctions
  #include half_norm_fcns.stanfunctions
}

data {
  int<lower=1> N;
  vector<lower=0>[N] fwi;
}

transformed data {
  int<lower=1> K = 2;
}

parameters {
  simplex[K] theta; // mixing proportions
  
  // egpd params
  real xi_prime;
  real kappa_prime;
  real sigma_prime;
  
  // // normal distribution params
  // real<lower=0> mu;
  // real<lower=0> sigma_norm;
  // exponential distribution params
  real<lower=0> scale;
}

transformed parameters{
  real<lower=-1> xi = log1p_exp(xi_prime) - 1;
  // real<lower=0> xi = log1p_exp(xi_prime);
  real<lower=0> kappa = exp(kappa_prime);
  real<lower=0> sigma_egpd;
  
  if(xi < 0) {
    real min_sigma_fwi = max(fwi) * (-xi);
    sigma_egpd = min_sigma_fwi + exp(sigma_prime);
  } else {
    sigma_egpd = exp(sigma_prime);
  }
  
}

model {
  vector[K] log_theta = log(theta);  // cache log calculation
  
  xi_prime ~ std_normal();
  kappa_prime ~ std_normal();
  sigma_prime ~ normal(0, 4);
  // sigma_egpd ~ exponential(1);
  
  // mu ~ normal(9, 10);
  // sigma_norm ~ lognormal(0, 2);
  scale ~ lognormal(0,2);
  
  for (n in 1:N) {
    vector[K] lps = log_theta;
    lps[1] += egpd_lpdf(fwi[n] | sigma_egpd, xi, kappa);
    // lps[2] += half_norm_lpdf(fwi[n] | sigma_norm);
    // lps[2] += normal_lpdf(fwi[n] | mu, sigma_norm);
    lps[2] += exponential_lpdf(fwi[n] | scale);
    target += log_sum_exp(lps);
  }
}
