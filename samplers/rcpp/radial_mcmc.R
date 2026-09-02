# =============================================================================
# Fit the radial component of the gauge-dependence model across all 200
# simulation-study datasets, using the adaptive Metropolis-Hastings sampler
# from the gaugeDependence package. The radial model conditions on the angles
# W and the gauge threshold r0(W), under either a truncated or censored
# likelihood.
#
# Called by: shell_scripts/run_radial_mcmc.sh
# Inputs:    data/{dep_type}/{dep_level}_{data_num}.json
# Outputs:   samplers/rcpp/radial_mcmc_fits/{dep_type}/
#              {gauge}_{likelihood}_{dep_level}_{data_num}.qs
#
# Command-line args:
#   1. dep_type   -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level  -- dependence strength ("low", "mid", "high")
#   3. gauge      -- gauge function ("gauss", "logistic", "inv_log", "asym_log",
#                                    "dirichlet", "rectangular")
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

# Starting values: the Dirichlet gauge has two shape parameters, so it needs a
# length-3 vector (alpha + two shapes); all other gauges need alpha + one dep.
if (gauge != "dirichlet") {
  starting_vals <- c(rgamma(1, 4, 2), runif(1))
} else {
  starting_vals <- c(rgamma(1, 4, 2), abs(rt(1, 4, ncp = 0)) * 4,
                                      abs(rt(1, 4, ncp = 0)) * 2)
}

# Fit each dataset in turn, saving one MCMC object per dataset.
for (data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  data <- RcppSimdJson::fload(datafile)

  w   <- data$W      # angular component
  r   <- data$R      # radial component
  r0w <- data$r0_w   # gauge threshold evaluated at each angle
  idx <- data$idx    # indices of the observations above threshold

  # The truncated likelihood is defined only over exceedances, so subset to
  # the above-threshold observations; the censored likelihood keeps them all.
  if (likelihood == "trunc") {
    w   <- w[idx]
    r   <- r[idx]
    r0w <- r0w[idx]
  }

  results <- radial_adaptive_mh(
    radii           = r,
    r0w             = r0w,
    angles          = w,
    starting_theta  = starting_vals,
    likelihood_type = likelihood,
    gauge_type      = gauge,
    n_updates       = 15000,
    update_freq     = 250,
    n_burnin        = 5000,
    n_thin          = 5,
    adapt_cov       = TRUE
  )

  qsave(x = results,
        file = sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs",
                       dep_type, gauge, likelihood, dep_level, data_num))
  print(paste0("Successfully saved MCMC fit for dataset number: ", data_num))
}
