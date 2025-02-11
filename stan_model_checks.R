library(cmdstanr)
model <- cmdstan_model("samplers/stan/marg_transform/fwi_transform_mix.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)
