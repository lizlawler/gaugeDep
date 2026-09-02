# =============================================================================
# Runs NIMBLE MCMC for the angular mixture model (Dirichlet process mixture
# of Betas, parameterised via a stick-breaking prior) across a specified
# range of simulation study datasets. The dataset range is controlled by
# command-line arguments so this script can be parallelised across jobs
# on an HPC cluster.
#
# Called by: shell_scripts/run_angle_mix_mcmc.sh
# Inputs:    data/{dep_type}/{dep_level}_{i}.json
# Outputs:   samplers/nimble/ang_mix_mcmc_fits/{dep_type}/{dep_level}_{i}.qs
#
# Command-line args:
#   1. dep_type   -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level  -- dependence strength ("low", "mid", "high")
#   3. data_start -- first dataset index to process
#   4. data_end   -- last dataset index to process
#
# Model overview:
#   Stick-breaking prior on mixture weights (GEM / DP prior):
#     v_i ~ Beta(1, alpha),  alpha ~ Gamma(2, 2)
#     probs <- stick_breaking(v)
#   Component parameters (mean/precision parameterisation of Beta):
#     beta_tau ~ Exp(1/8)          [precision hyperprior]
#     mu_k ~ Uniform(0, 1)         [component mean]
#     tau_k ~ InvGamma(2, beta_tau)[component precision]
#     alpha_k = mu_k * tau_k,  beta_k = (1 - mu_k) * tau_k
#   Likelihood:
#     z_i ~ Categorical(probs),  w_i ~ Beta(alpha_{z_i}, beta_{z_i})
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type   <- args[1]
dep_level  <- args[2]
data_start <- as.integer(args[3])
data_end   <- as.integer(args[4])

library(nimble)
library(tidyr)
library(dplyr)
library(qs)

# Stick-breaking Dirichlet process mixture of Betas (L = truncation level)
sb_code <- nimbleCode({
  for (i in 1:(L - 1)) {
    v[i] ~ dbeta(1, alpha)
  }
  alpha    ~ dgamma(2, 2)
  beta_tau ~ dexp(1/8)
  probs[1:L] <- stick_breaking(v[1:(L - 1)])

  for (i in 1:L) {
    mustar[i]    ~ dunif(0, 1)
    taustar[i]   ~ dinvgamma(2, beta_tau)
    alphastar[i] <- mustar[i] * taustar[i]
    betastar[i]  <- (1 - mustar[i]) * taustar[i]
  }
  for (i in 1:N) {
    z[i]  ~ dcat(probs[1:L])
    w[i]  ~ dbeta(alphastar[z[i]], betastar[z[i]])
  }
})

for (i in data_start:data_end) {
  datafile <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, i)
  data     <- RcppSimdJson::fload(datafile)

  sb_data      <- list(w = data$W)
  sb_constants <- list(N = length(sb_data$w), L = 10)

  # Random initialisations for mixture parameters
  sb_inits <- list(
    mustar   = runif(sb_constants$L, 0, 1),
    taustar  = rinvgamma(sb_constants$L, 1, 1),
    z        = sample(1:5, size = sb_constants$N, replace = TRUE),
    v        = rbeta(sb_constants$L, 1, 1),
    alpha    = 1,
    beta_tau = 8
  )

  # Build, compile, configure, and run NIMBLE MCMC
  sb_model   <- nimbleModel(sb_code, sb_constants, sb_data, sb_inits)
  cmodel     <- compileNimble(sb_model, resetFunctions = TRUE)
  conf_model <- configureMCMC(sb_model)

  # Monitor mixture parameters needed for downstream loglik computation
  conf_model$addMonitors("alphastar", "betastar", "taustar", "probs", "z")
  modelMCMC  <- buildMCMC(conf_model)
  cmodelMCMC <- compileNimble(modelMCMC, project = sb_model, resetFunctions = TRUE)

  results <- runMCMC(cmodelMCMC,
                     niter   = 15000,
                     nburnin = 5000,
                     thin    = 5,
                     inits   = sb_inits,
                     nchains = 1)

  qsave(x = results,
        file = sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
                       dep_type, dep_level, i))
  print(paste0("Successfully saved MCMC stick-breaking fit for dataset number: ", i))
}
