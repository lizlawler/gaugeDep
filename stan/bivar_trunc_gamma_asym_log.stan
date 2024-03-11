functions {
  #include gauge_fcns.stanfunctions
  #include truncGamma.stanfunctions
}
#include data_params.stan
transformed parameters {
  vector[N] beta;
  for (n in 1:N) {
    beta[n] = asym_log_gauge(W[n], dep);
  }
}
#include model_gq.stan
