#include fcns_data_params_trunc.stan
transformed parameters {
  // vector[n0] beta;
  vector[n0_ctau] beta;
  // for (n in 1:n0) {
  for (n in 1:n0_ctau) {
    beta[n] = gauss_gauge(W_trunc[n], dep);
  }
}
#include model_gq_trunc.stan
