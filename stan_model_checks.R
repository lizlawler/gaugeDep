library(cmdstanr)
model <- cmdstan_model("stan/radial_angular/bivar_cens_marg_gauss_mix_betas.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)
