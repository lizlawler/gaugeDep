# =============================================================================
# Runs the BezELS (Bezier Extreme Level Set) MCMC sampler for the radial
# component across all 200 simulation study datasets. BezELS provides a
# nonparametric alternative to the parametric gauge-function radial model and
# serves as a competitor method in the prediction task.
#
# Called by: shell_scripts/run_bezels_mcmc.sh
# Inputs:    data/{dep_type}/{dep_level}_{data_num}.json
# Outputs:   samplers/bezels/radial_bezels_fits/{dep_type}/{dep_level}_{data_num}.qs
#
# Command-line args:
#   1. dep_type  -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level -- dependence strength ("low", "mid", "high")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type  <- args[1]
dep_level <- args[2]

library(dplyr)
library(tidyr)
library(qs)
library(BezELS)

for (data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  data     <- RcppSimdJson::fload(datafile)
  idx      <- data$idx

  # Restrict to threshold exceedances (BezELS is fit only above the threshold)
  w   <- data$W[idx]
  r   <- data$R[idx]
  r0w <- data$r0_w[idx]

  # N is the total sample size (used to recover the unconditional scale).
  results <- fit_mcmc_bezier(
    N           = data$n0,
    r           = r,
    w           = w,
    r_0         = r0w,
    iters       = 11000,
    burn        = 100,
    traceplot   = FALSE,
    print.every = 100
  )

  qsave(x = results,
        file = sprintf("samplers/bezels/radial_bezels_fits/%s/%s_%s.qs",
                       dep_type, dep_level, data_num))
  print(paste0("Successfully saved BezELS fit for dataset number: ", data_num))
}
