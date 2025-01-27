library(cmdstanr)
model <- cmdstan_model("samplers/stan/marg_transform/fire_transform.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)
