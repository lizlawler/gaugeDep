#include fcns_data_params_cens.stan
model {
  alpha ~ gamma(4, 2);
  dep ~ uniform(0, 1);
  
  for (n in 1:N) {
    target += cens_gamma_lpdf(R[n] | r0_w[n], alpha, asym_log_gauge(W[n], dep));
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = cens_gamma_lpdf(R[n] | r0_w[n], alpha, asym_log_gauge(W[n], dep));
  }
}

