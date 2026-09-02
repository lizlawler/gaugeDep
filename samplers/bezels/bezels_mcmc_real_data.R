# =============================================================================
# Runs the BezELS MCMC sampler on real fire weather data (one station at a
# time). Currently hardcoded to "redstone"; update data_type to switch stations.
#
# Run interactively (no command-line args).
# Inputs:    data/raw/{data_type}_expo.qs
# Outputs:   samplers/bezels/radial_bezels_fits/real_data/{data_type}.qs
# =============================================================================

library(dplyr)
library(tidyr)
library(qs)
library(BezELS)

data_type <- "redstone"  # update to switch stations ("friendmtn", "redstone")

data <- qread(sprintf("data/raw/%s_expo.qs", data_type))
idx  <- data$idx

# Restrict to threshold exceedances
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
      file = sprintf("samplers/bezels/radial_bezels_fits/real_data/%s.qs", data_type))
print(paste0("Successfully saved BezELS fit on fire data for ", data_type))
