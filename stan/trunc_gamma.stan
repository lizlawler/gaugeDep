/*functions {
  vector trunc_gamma_rng(int n, real L, real alpha, real beta) {
    // returns a vector length n
    real cst = gamma_cdf(L | alpha, beta);
    vector[n] rng_val;
    vector[n] a = rep_vector(0, n);
    vector[n] b = rep_vector(1, n);
    array[n] real u = uniform_rng(a, b);
    for (i in 1:n) {
      real u_adj = u[i] * (1-cst);
      rng_val[i] = gamma_icdf(u_adj, alpha, beta);
    }
    return rng_val;
  }
}
*/
data {
  int<lower=1> N;
  vector[N] x;
  real L;
}
parameters {
  /*real<lower = 0, upper=min(x)> L;*/
  real<lower=0> alpha;
  real<lower=0> beta;
}
model {
  /*L ~ exponential(0.1);*/
  alpha ~ exponential(0.01);
  beta ~ exponential(0.01);
  x ~ gamma(alpha, beta) T[L, ];
}
/*generated quantities {
  vector[N] pred_x;
  pred_x = trunc_gamma_rng(N, L, alpha, beta);
  gamma_q()
}*/
