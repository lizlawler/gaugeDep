#include fcns_data_params_cens.stan
transformed parameters {
  vector[N] beta;
  for (n in 1:N) {
    beta[n] = asym_log_gauge(W[n], dep);
  }
}
#include model_gq_cens.stan
