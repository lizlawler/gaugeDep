// shared functions, data input values and shared parameter declaration
functions {
  #include gauge_fcns.stanfunctions
  #include likelihoodGamma.stanfunctions
}

data {
  int<lower=1> N;
  array[N] real<lower=0> R;
  array[N] real<lower=0, upper=1> W;
  int<lower=1> n0;
  array[N] real<lower=0> r0_w;
  array[n0] int<lower=1> idx;
  // int<lower=1> n0_ctau;
  // array[N] real<lower=0> r0_w_ctau;
  // array[n0_ctau] int<lower=1> idx_ctau;
}

transformed data {
  array[n0] real<lower=0> R_trunc = R[idx];
  array[n0] real<lower=0, upper=1> W_trunc = W[idx];
  array[n0] real<lower=0> r0_w_trunc = r0_w[idx];
  // array[n0_ctau] real<lower=0> R_trunc = R[idx_ctau];
  // array[n0_ctau] real<lower=0, upper=1> W_trunc = W[idx_ctau];
  // array[n0_ctau] real<lower=0> r0_w_trunc = r0_w_ctau[idx_ctau];
}

parameters {
  real<lower=0> alpha;
  real<lower=0, upper =1> dep;
}