args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
data_batch <- args[3]

library(nimble)
library(tidyr)
library(dplyr)
library(qs)

sb_code <- nimbleCode({
  for(i in 1:(L-1)){
    v[i] ~ dbeta(1, alpha)
  }
  alpha ~ dgamma(2, 2)
  beta_tau ~ dexp(1/8)
  probs[1:L] <- stick_breaking(v[1:(L-1)])
  for(i in 1:L) {
    mustar[i] ~ dunif(0, 1)
    taustar[i] ~ dinvgamma(2, beta_tau)
    alphastar[i] <- mustar[i] * taustar[i]
    betastar[i] <- (1-mustar[i]) * taustar[i]
  }
  for(i in 1:N) {
    z[i] ~ dcat(probs[1:L])
    w[i] ~ dbeta(alphastar[z[i]], betastar[z[i]])
  }
})

data_basename <- paste0("data/", dep_type, "/", dep_level, "_")
batch_nums <- if(data_batch == 1) {
  1:25
} else if(data_batch == 2) {
  26:50
} else if(data_batch == 3) {
  51:75
} else {
  76:100
}

for(i in seq_along(batch_nums)) {
  data_num <- batch_nums[i]
  datafile <- paste0(data_basename, data_num, ".json")
  sb_data <- list(w = RcppSimdJson::fload(datafile)$W)
  sb_constants <- list(N = length(sb_data$w), L=10)
  sb_inits <- list(mustar = runif(sb_constants$L, 0, 1), 
                   taustar = rinvgamma(sb_constants$L, 1, 1), 
                   z = sample(1:10, size = sb_constants$N, replace = TRUE),
                   v  = rbeta(sb_constants$L, 1, 1),
                   alpha = 1,
                   beta_tau = 8)
  
  sb_model <- nimbleModel(sb_code, sb_constants, sb_data, sb_inits)
  cmodel <- compileNimble(sb_model, resetFunctions = TRUE)
  conf_model <- configureMCMC(sb_model)
  conf_model$addMonitors("alphastar", "betastar", "taustar", "probs", "z")
  modelMCMC <- buildMCMC(conf_model)
  
  cmodelMCMC <- compileNimble(modelMCMC, project = sb_model, resetFunctions = TRUE)
  results <- runMCMC(cmodelMCMC, 
                     niter = 15000, nburnin = 5000, thin = 5, 
                     inits = sb_inits, nchains = 1)

  qsave(x = results, file = paste0("samplers/nimble/sb_mcmc_fits/", dep_type, "/", 
                                   dep_level, "_", data_num, ".qs"))
  print(paste0("Successfully saved MCMC stick breaking fit for dataset number: ", data_num))
}
