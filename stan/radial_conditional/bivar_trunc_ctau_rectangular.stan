#include fcns_data_params_trunc_ctau.stan
model {
  alpha ~ gamma(4, 2);
  dep ~ uniform(0, 1);
  
  for (n in 1:n0_ctau) {
    target += trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, rectangular_gauge(W_trunc[n], dep));
  }
}

generated quantities {
  vector[n0_ctau] log_lik;
  for (n in 1:n0_ctau) {
    log_lik[n] = trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, rectangular_gauge(W_trunc[n], dep));
  }
}
