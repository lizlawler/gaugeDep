# =============================================================================
# Compute the posterior pointwise log-likelihood of the star-shaped angular
# density on real fire-weather data. Output is an (n_iter x n_obs) matrix
# consumed by the loo package to form stacking / pseudo-BMA weights.
#
# Called by: shell_scripts/local_machine/run_angular.sh (or similar)
# Inputs:    data/raw/{data_type}_expo.qs
#            samplers/rcpp/ang_star_mcmc_fits/real_data/{data_type}_{gauge}_2chains.qs
# Outputs:   samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/
#              {data_type}_{gauge}.qs
#
# Command-line args:
#   1. data_type -- station identifier ("friendmtn" or "redstone")
#   2. gauge     -- gauge function ("gauss", "logistic", "inv_log", etc.)
#
# Note: to run interactively, set data_type and gauge manually.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_type <- args[1]
gauge     <- args[2]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

# The multi-chain object stores combined draws in posterior::draws_df format.
params <- qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/%s_%s_2chains.qs",
                        data_type, gauge))
data <- qread(sprintf("data/raw/%s_expo.qs", data_type))
w    <- data$W

results <- angular_loglik(
  angles           = w,
  dim              = 2,
  posterior_params = params,
  gauge_type       = gauge
)

qsave(x = results,
      file = sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/%s_%s.qs",
                     data_type, gauge))
print(sprintf(
  "Successfully saved posterior pointwise loglikelihood of star-shaped density for %s, %s gauge",
  data_type, gauge
))
