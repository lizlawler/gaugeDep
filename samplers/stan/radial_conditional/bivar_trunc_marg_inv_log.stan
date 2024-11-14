#include fcns_data_params_trunc_marg.stan
model {
  alpha ~ gamma(4, 2);
  dep ~ uniform(0, 1);
  
  for (n in 1:n0) {
    target += trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, inv_log_gauge(W_trunc[n], dep));
  }
}

generated quantities {
  vector[n0] log_lik;
  for (n in 1:n0) {
    log_lik[n] = trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, inv_log_gauge(W_trunc[n], dep));
  }
}
