args <- commandArgs(trailingOnly=TRUE)
gauge <- args[1]
likelihood <- args[2]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

data_type <- "redstone" 

data <- qread(sprintf("data/%s_expo.qs", data_type))
params <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/real_data/%s_%s_%s.qs", data_type, gauge, likelihood))$samples
params <- params[,1:(ncol(params) - 1)]
w <- data$W
r <- data$R
r0w <- data$r0_w
idx <- data$idx
if(likelihood == "trunc") {
  w <- w[idx]
  r <- r[idx]
  r0w <- r0w[idx]
}

results <- radial_loglik(radii = r, threshold = r0w, angles = w,
                         posterior_params = params,
                         likelihood_type = likelihood,
                         gauge_type = gauge)
qsave(x = results, file = sprintf("samplers/rcpp/radial_mcmc_fits/real_data/pw_loglik/%s_%s_%s.qs", data_type, gauge, likelihood))
print(sprintf("Successfully saved posterior pointwise loglikelihood for radial density of %s, with %s gauge and %s likelihood", data_type, gauge, likelihood))
