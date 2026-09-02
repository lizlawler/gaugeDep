# =============================================================================
# Compute the posterior pointwise log-likelihood of the radial component for
# every simulation-study dataset. Each output is an (n_iter x n_obs) matrix,
# consumed downstream by the loo package to form stacking / pseudo-BMA weights.
#
# Called by: shell_scripts/run_radial_loglik_calc.sh
# Inputs:    data/{dep_type}/{dep_level}_{data_num}.json
#            samplers/rcpp/radial_mcmc_fits/{dep_type}/
#              {gauge}_{likelihood}_{dep_level}_{data_num}.qs
# Outputs:   samplers/rcpp/radial_mcmc_fits/{dep_type}/pw_loglik/
#              {gauge}_{likelihood}_{dep_level}_{data_num}.qs
#
# Command-line args:
#   1. dep_type   -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level  -- dependence strength ("low", "mid", "high")
#   3. gauge      -- gauge function ("gauss", "logistic", "inv_log", etc.)
#   4. likelihood -- likelihood type ("trunc" or "cens")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type   <- args[1]
dep_level  <- args[2]
gauge      <- args[3]
likelihood <- args[4]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

for (data_num in 1:200) {
  datafile   <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  paramsfile <- sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs",
                        dep_type, gauge, likelihood, dep_level, data_num)

  data <- RcppSimdJson::fload(datafile)

  # Drop the final column of the sample matrix, which holds the sampler's
  # log-posterior bookkeeping rather than a model parameter.
  params <- qread(paramsfile)$samples
  params <- params[, 1:(ncol(params) - 1)]

  w   <- data$W
  r   <- data$R
  r0w <- data$r0_w
  idx <- data$idx

  # Match the subsetting used when the model was fit: truncated uses
  # exceedances only, censored uses all observations.
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
        file = sprintf("samplers/rcpp/radial_mcmc_fits/%s/pw_loglik/%s_%s_%s_%s.qs",
                       dep_type, gauge, likelihood, dep_level, data_num))
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}
