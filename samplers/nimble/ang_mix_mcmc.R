args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
data_start <- args[3]
data_end <- args[4]

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

for(i in data_start:data_end) {
  datafile <- sprintf("data/%s/%s_%s.json",
                      dep_type, dep_level, i)
  data <- RcppSimdJson::fload(datafile)
  sb_data <- list(w = data$W)
  sb_constants <- list(N = length(sb_data$w), L=10)
  sb_inits <- list(mustar = runif(sb_constants$L, 0, 1), 
                   taustar = rinvgamma(sb_constants$L, 1, 1), 
                   z = sample(1:5, size = sb_constants$N, replace = TRUE),
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

  qsave(x = results, file = sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
                                    dep_type, dep_level, i))
  print(paste0("Successfully saved MCMC stick breaking fit for dataset number: ", i))
}
