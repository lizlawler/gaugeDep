// shared functions, data input values and shared parameter declaration
functions {
  #include gauge_fcns.stanfunctions
  #include likelihoodGamma.stanfunctions
}

data {
  int<lower=1> N;
  array[N] real<lower=0> R;
  array[N] real<lower=0, upper=1> W;
  array[N] real<lower=0> r0_w;
  array[N] real<lower=0> r0_w_ctau;
}

parameters {
  real<lower=0> alpha;
  real<lower=0, upper =1> dep;
}