// same model and gen_quant declarations across differing gauge functions
model {
  alpha ~ gamma(4, 2);
  dep ~ uniform(0, 1);
  
  for (n in 1:n0) {
  // for (n in 1:n0_ctau) {
    target += trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, beta[n]);
  }
}

generated quantities {
  vector[n0] log_lik;
  // vector[n0_ctau] log_lik;
  for (n in 1:n0) {
  // for (n in 1:n0_ctau) {
    log_lik[n] = trunc_gamma_lpdf(R_trunc[n] | r0_w_trunc[n], alpha, beta[n]);
  }
}
