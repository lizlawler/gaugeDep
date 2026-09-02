# =============================================================================
# Computes the posterior pointwise log-likelihood for the angular mixture
# (Dirichlet process mixture of Betas) model across all simulation study
# datasets. Output is used downstream for stacking weight computation via
# the `loo` package.
#
# Called by: shell_scripts/run_angle_mix_loglik_calc.sh
# Inputs:    data/{dep_type}/{dep_level}_{data_num}.json
#            samplers/nimble/ang_mix_mcmc_fits/{dep_type}/{dep_level}_{data_num}.qs
# Outputs:   samplers/nimble/ang_mix_mcmc_fits/{dep_type}/pw_loglik/{dep_level}_{data_num}.qs
#
# Command-line args:
#   1. dep_type  -- dependence structure (e.g., "gauss", "logistic")
#   2. dep_level -- dependence strength (e.g., "low", "mid", "high")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type  <- args[1]
dep_level <- args[2]

library(nimble)
library(tidyr)
library(dplyr)
library(qs)

# Shared helper functions (mix_lpdf, angular_loglik)
source("samplers/nimble/angular_mix_helpers.R")

for (data_num in 1:200) {
  datafile   <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  paramsfile <- sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
                        dep_type, dep_level, data_num)

  data   <- RcppSimdJson::fload(datafile)
  params <- qread(paramsfile)
  w      <- data$W

  # Compute (n_iter x n_obs) matrix of pointwise log-likelihoods
  results <- angular_loglik(angles = w, posterior_params = params)

  savename <- sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/pw_loglik/%s_%s.qs",
                      dep_type, dep_level, data_num)
  qsave(x = results, file = savename)
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}
