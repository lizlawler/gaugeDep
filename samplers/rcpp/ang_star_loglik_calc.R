# =============================================================================
# Compute the posterior pointwise log-likelihood of the star-shaped angular
# density for every simulation-study dataset. Each output is an
# (n_iter x n_obs) matrix, consumed downstream by the loo package to form
# stacking / pseudo-BMA weights.
#
# Called by: shell_scripts/run_angle_star_loglik_calc.sh
# Inputs:    data/{dep_type}/{dep_level}_{data_num}.json
#            samplers/rcpp/ang_star_mcmc_fits/{dep_type}/{gauge}_{dep_level}_{data_num}.qs
# Outputs:   samplers/rcpp/ang_star_mcmc_fits/{dep_type}/pw_loglik/
#              {gauge}_{dep_level}_{data_num}.qs
#
# Command-line args:
#   1. dep_type  -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level -- dependence strength ("low", "mid", "high")
#   3. gauge     -- gauge function ("gauss", "logistic", "inv_log", etc.)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type  <- args[1]
dep_level <- args[2]
gauge     <- args[3]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

print(paste0("dep_type = ",  dep_type))
print(paste0("dep_level = ", dep_level))
print(paste0("gauge = ",     gauge))

for (data_num in 1:200) {
  datafile   <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  paramsfile <- sprintf("samplers/rcpp/ang_star_mcmc_fits/%s/%s_%s_%s.qs",
                        dep_type, gauge, dep_level, data_num)

  data   <- RcppSimdJson::fload(datafile)
  w      <- data$W

  # Drop the final column (sampler log-posterior bookkeeping, not a parameter).
  params <- qread(paramsfile)$samples
  params <- params[, 1:(ncol(params) - 1)]

  results <- angular_loglik(
    angles           = w,
    dim              = 2,
    posterior_params = params,
    gauge_type       = gauge
  )

  savename <- sprintf("samplers/rcpp/ang_star_mcmc_fits/%s/pw_loglik/%s_%s_%s.qs",
                      dep_type, gauge, dep_level, data_num)
  qsave(x = results, file = savename)
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}
