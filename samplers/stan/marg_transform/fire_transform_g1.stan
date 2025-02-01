functions {
  #include gpd_fcns.stanfunctions
  #include g1_fcns.stanfunctions
}

data {
  int<lower=1> N;
  int<lower=1> D;
  vector<lower=0>[N] erc;
  vector<lower=0>[N] fwi;
}

parameters {
  vector[D] xi_prime;
  vector[D] kappa_prime;
  vector[D] sigma;
}

transformed parameters{
  vector<lower=-1>[D] xi = log1p_exp(xi_prime) - 1;
  vector<lower=0>[D] kappa = exp(kappa_prime);
  // vector<lower=0>[D] sigma;
  // // vector<lower=0>[D] sigma = exp(sigma_prime);
  // 
  // real max_sigma_erc = max(erc) * (-xi[1]);  
  // real max_sigma_fwi = max(fwi) * (-xi[2]);  
  // 
  // sigma[1] = max_sigma_erc * exp(sigma_prime[1]);
  // sigma[2] = max_sigma_fwi * exp(sigma_prime[2]);
}

model {
  to_vector(xi_prime) ~ std_normal();
  to_vector(kappa_prime) ~ std_normal();
  to_vector(sigma) ~ exponential();
  // likelihood
  for(n in 1:N) {
    target += egpd_lpdf(erc[n] | sigma[1], xi[1], kappa[1]);
    target += egpd_lpdf(fwi[n] | sigma[2], xi[2], kappa[2]);
  }
}
