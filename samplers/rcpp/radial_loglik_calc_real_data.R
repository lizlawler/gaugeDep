# =============================================================================
# Compute the posterior pointwise log-likelihood of the radial component on
# real fire-weather data. Output is an (n_iter x n_obs) matrix consumed by the
# loo package to form stacking / pseudo-BMA weights.
#
# Called by: shell_scripts/local_machine/run_joint.sh (or similar)
# Inputs:    data/raw/{data_type}_expo.qs
#            samplers/rcpp/radial_mcmc_fits/real_data/
#              {data_type}_{gauge}_{likelihood}_2chains.qs
# Outputs:   samplers/rcpp/radial_mcmc_fits/real_data/pw_loglik/
#              {data_type}_{gauge}_{likelihood}.qs
#
# Command-line args:
#   1. data_type  -- station identifier ("friendmtn" or "redstone")
#   2. gauge      -- gauge function ("gauss", "logistic", "inv_log", etc.)
#   3. likelihood -- likelihood type ("trunc" or "cens")
#
# Note: to run interactively, set data_type, gauge, and likelihood manually.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_type  <- args[1]
gauge      <- args[2]
likelihood <- args[3]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

data <- qread(sprintf("data/raw/%s_expo.qs", data_type))

# The multi-chain object is already in draws_df format (posterior package), so
# unlike the simulation-study scripts there is no bookkeeping column to drop.
params <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/real_data/%s_%s_%s_2chains.qs",
                        data_type, gauge, likelihood))

w   <- data$W
r   <- data$R
r0w <- data$r0_w
idx <- data$idx

if (likelihood == "trunc") {
  w   <- w[idx]
  r   <- r[idx]
  r0w <- r0w[idx]
}

results <- radial_loglik(
  radii            = r,
  threshold        = r0w,
  angles           = w,
  posterior_params = params,
  likelihood_type  = likelihood,
  gauge_type       = gauge
)

qsave(x = results,
      file = sprintf("samplers/rcpp/radial_mcmc_fits/real_data/pw_loglik/%s_%s_%s.qs",
                     data_type, gauge, likelihood))
print(sprintf(
  "Successfully saved posterior pointwise loglikelihood for radial density of %s, %s gauge, %s likelihood",
  data_type, gauge, likelihood
))
