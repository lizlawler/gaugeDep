# =============================================================================
# Computes the posterior pointwise log-likelihood for the angular mixture
# (Dirichlet process mixture of Betas) model on real fire weather data.
# Output is used downstream for stacking weight computation via the `loo`
# package.
#
# Called by: shell_scripts/local_machine/run_angular.sh (or similar)
# Inputs:    data/raw/{data_type}_expo.qs
#            samplers/nimble/ang_mix_mcmc_fits/real_data/{data_type}_2chains.qs
# Outputs:   samplers/nimble/ang_mix_mcmc_fits/real_data/pw_loglik/{data_type}.qs
#
# Command-line args:
#   1. data_type -- station identifier (e.g., "friendmtn", "redstone")
#
# Note: to run interactively for a single station, set data_type manually
#       (e.g., data_type <- "friendmtn") rather than via command-line args.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_type <- args[1]

library(nimble)
library(tidyr)
library(dplyr)
library(qs)

# Shared helper functions (mix_lpdf, angular_loglik)
source("samplers/nimble/angular_mix_helpers.R")

data   <- qread(sprintf("data/raw/%s_expo.qs", data_type))
params <- qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/%s_2chains.qs", data_type))
w      <- data$W

# Compute (n_iter x n_obs) matrix of pointwise log-likelihoods
results <- angular_loglik(angles = w, posterior_params = params)

qsave(x = results,
      file = sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/pw_loglik/%s.qs", data_type))
print(sprintf("Successfully saved posterior pointwise loglikelihood for angular mixture density of %s",
              data_type))
