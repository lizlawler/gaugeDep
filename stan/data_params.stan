// shared data input values and shared parameter declaration
data {
  int<lower=0> N;
  vector[N] R;
  vector[N] W;
  vector[N] r0_w;
}

parameters {
  real<lower=0> alpha;
  real<lower=0, upper =1> dep;
}