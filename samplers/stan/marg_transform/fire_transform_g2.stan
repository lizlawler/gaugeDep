functions {
  #include gpd_fcns.stanfunctions
  #include g2_fcns.stanfunctions
}

data {
  int<lower=1> N;
  int<lower=1> D;
  vector<lower=0>[N] erc;
  vector<lower=0>[N] fwi;
}

parameters {
  vector[D] log_xi;
  vector[D] log_kappa1;
  vector[D] log_kappa2;
  vector[D] log_sigma;
  vector<lower=0,upper=1>[D] probs;
}

transformed parameters{
  vector<lower=0>[D] xi = exp(log_xi);
  vector<lower=0>[D] kappa1 = exp(log_kappa1);
  vector<lower=0>[D] kappa2 = exp(log_kappa2);
  vector<lower=0>[D] sigma = exp(log_sigma);
}

model {
  to_vector(probs) ~ uniform(0, 1);
  to_vector(log_xi) ~ std_normal();
  to_vector(log_kappa1) ~ std_normal();
  to_vector(log_kappa2) ~ std_normal();
  to_vector(log_sigma) ~ std_normal();
  // likelihood
  for(n in 1:N) {
    target += egpd_lpdf(erc[n] | sigma[1], xi[1], kappa1[1], kappa2[1], probs[1]);
    target += egpd_lpdf(fwi[n] | sigma[2], xi[2], kappa1[2], kappa2[2], probs[2]);
  }
}
