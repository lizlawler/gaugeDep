args <- commandArgs(trailingOnly=TRUE)
gauge <- args[1]
likelihood <- args[2]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

data_type <- "redstone"

if(gauge != "dirichlet") {
  starting_vals <- c(rgamma(1, 4, 2), runif(1))
} else {
  starting_vals <- c(rgamma(1, 4, 2), abs(rt(1, 4,ncp = 0))*4, abs(rt(1, 4,ncp = 0))*2)
}

data <- qread(sprintf("data/%s_expo.qs", data_type))
w <- data$W
r <- data$R
r0w <- data$r0_w
idx <- data$idx
if(likelihood == "trunc") {
  w <- w[idx]
  r <- r[idx]
  r0w <- r0w[idx]
}

results <- radial_adaptive_mh(radii = r, r0w = r0w, angles = w,
                              starting_theta = starting_vals,
                              likelihood_type = likelihood,
                              gauge_type = gauge,
                              n_updates = 15000, 
                              update_freq = 250, 
                              n_burnin = 5000,
                              n_thin = 5,
                              adapt_cov = TRUE)
qsave(x = results, file = sprintf("samplers/rcpp/radial_mcmc_fits/real_data/%s_%s_%s.qs", data_type, gauge, likelihood))
print(sprintf("Successfully saved radial MCMC fit on fire data for %s, with %s gauge and %s likelihood", data_type, gauge, likelihood))


