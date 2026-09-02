# =============================================================================
# Fit the star-shaped angular density model across all 200 simulation-study
# datasets, using the adaptive Metropolis-Hastings sampler from the
# gaugeDependence package. The angular model is fit independently of the radial
# component and conditions only on the observed angles W.
#
# Called by: shell_scripts/run_angle_star_mcmc.sh
# Inputs:    data/{dep_type}/{dep_level}_{data_num}.json
# Outputs:   samplers/rcpp/ang_star_mcmc_fits/{dep_type}/
#              {gauge}_{dep_level}_{data_num}.qs
#
# Command-line args:
#   1. dep_type  -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level -- dependence strength ("low", "mid", "high")
#   3. gauge     -- gauge function ("gauss", "logistic", "inv_log", "asym_log",
#                                   "dirichlet", "rectangular")
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

# Starting values: the Dirichlet gauge has two shape parameters; every other
# gauge has a single scalar dependence parameter.
if (gauge != "dirichlet") {
  starting_vals <- runif(1)
} else {
  starting_vals <- c(abs(rt(1, 4, ncp = 0)) * 4,
                     abs(rt(1, 4, ncp = 0)) * 2)
}

# Fit each dataset in turn. Only the angles W are needed here; the radial
# component plays no part in the angular model.
for (data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  data <- RcppSimdJson::fload(datafile)
  w <- data$W

  results <- angular_mcmc(
    angles         = w,
    dim            = 2,
    starting_theta = starting_vals,
    gauge_type     = gauge,
    n_updates      = 15000,
    update_freq    = 250,
    n_burnin       = 5000,
    n_thin         = 5,
    adapt_cov      = TRUE
  )

  qsave(x = results,
        file = sprintf("samplers/rcpp/ang_star_mcmc_fits/%s/%s_%s_%s.qs",
                       dep_type, gauge, dep_level, data_num))
  print(paste0("Successfully saved MCMC fit for dataset number: ", data_num))
}
