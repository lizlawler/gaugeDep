// same model and gen_quant declarations across differing gauge functions
model {
  alpha ~ gamma(2, 4);
  dep ~ beta(1.5, 1.5);
  
  for (n in 1:N) {
    target += trunc_gamma_lpdf(R[n] | r0_w[n], alpha, beta[n]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = trunc_gamma_lpdf(R[n] | r0_w[n], alpha, beta[n]);
  }
}