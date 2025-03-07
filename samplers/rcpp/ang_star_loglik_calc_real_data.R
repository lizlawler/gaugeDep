args <- commandArgs(trailingOnly=TRUE)
gauge <- args[1]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

data_type <- "redstone"

paramsfile <- sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/%s_%s.qs", data_type, gauge)
data <- qread(sprintf("data/%s_expo.qs", data_type))
params <- qread(paramsfile)$samples
params <- params[,1:(ncol(params) - 1)]
w <- data$W

results <- angular_loglik(angles = w, dim = 2,
                          posterior_params = params, 
                          gauge_type = gauge)

qsave(x = results, file = sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/%s_%s.qs", data_type, gauge))
print(sprintf("Successfully saved posterior pointwise loglikelihood of star-shaped density for %s, %s gauge", data_type, gauge))


