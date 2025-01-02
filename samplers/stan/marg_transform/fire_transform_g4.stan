functions {
  #include gpd_fcns.stanfunctions
  #include g4_fcns.stanfunctions
}

data {
  int<lower=1> N;
  int<lower=1> D;
  vector<lower=0>[N] erc;
  vector<lower=0>[N] fwi;
}

parameters {
  vector[D] log_xi;
  vector[D] log_delta;
  vector[D] log_kappa;
  vector[D] log_sigma;
}

transformed parameters{
  vector<lower=0>[D] xi = exp(log_xi);
  vector<lower=0>[D] delta = exp(log_delta);
  vector<lower=0>[D] kappa = exp(log_kappa);
  vector<lower=0>[D] sigma = exp(log_sigma);
}

model {
  to_vector(log_xi) ~ std_normal();
  to_vector(log_delta) ~ std_normal();
  to_vector(log_kappa) ~ std_normal();
  to_vector(log_sigma) ~ std_normal();
  // likelihood
  for(n in 1:N) {
    target += egpd_lpdf(erc[n] | sigma[1], xi[1], delta[1], kappa[1]);
    target += egpd_lpdf(fwi[n] | sigma[2], xi[2], delta[2], kappa[2]);
  }
}
