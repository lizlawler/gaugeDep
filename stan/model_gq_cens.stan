// same model and gen_quant declarations across differing gauge functions
model {
  alpha ~ gamma(4, 2);
  dep ~ uniform(0, 1);
  
  for (n in 1:N) {
    target += cens_gamma_lpdf(R[n] | r0_w[n], alpha, beta[n]);
    // target += cens_gamma_lpdf(R[n] | r0_w_ctau[n], alpha, beta[n]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = cens_gamma_lpdf(R[n] | r0_w[n], alpha, beta[n]);
    // log_lik[n] = cens_gamma_lpdf(R[n] | r0_w_ctau[n], alpha, beta[n]);
  }
}
