library(cmdstanr)
model <- cmdstan_model("samplers/stan/marg_transform/fwi_transform_mix.stan", compile = TRUE)
model$check_syntax(pedantic = TRUE)
model$expose_functions()
curve(model$functions$trunc_expo_lpdf(x, 60, 0.05))

x <- seq(0, 100, length.out =1000)
lpdf <- rep(NA, length(x))
for(i in seq_along(x)) {
  lpdf[i] <- model$functions$trunc_expo_lpdf(x[i], 60, 0.05)
}
plot(x, exp(lpdf))
model$functions$trunc_expo_lpdf(75, 60, 0.05)

trunc_exp_r <- function(x, xmax, rate) {
  num <- dexp(x, rate = rate, log = TRUE)
  den <- pexp(xmax, rate = rate, log = TRUE)
  return(num - den)
}

exp(trunc_exp_r(75, 60, 0.05))

exp(model$functions$egpd_trunc_lpdf(1.5,1,0.5,0.75,2))
exp(model$functions$egpd_lpdf(1.5,0.5,0.75,2))
curve(dexp(x, 0.05), add = TRUE)
