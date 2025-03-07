args <- commandArgs(trailingOnly=TRUE)
gauge <- args[1]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

data_type <- "redstone"

if(gauge != "dirichlet") {
  starting_vals <- runif(1)
} else {
  starting_vals <- c(abs(rt(1, 4,ncp = 0))*4, abs(rt(1, 4,ncp = 0))*2)
}

data <- qread(sprintf("data/%s_expo.qs", data_type))
w <- data$W

results <- angular_mcmc(angles = w, dim = 2, 
                        starting_theta = starting_vals, 
                        gauge_type = gauge, 
                        n_updates = 15000, 
                        update_freq = 250, 
                        n_burnin = 5000, 
                        n_thin = 5, 
                        adapt_cov = TRUE)
qsave(x = results, file = sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/%s_%s.qs", data_type, gauge))
print(sprintf("Successfully saved star-shaped MCMC fit on % for %s gauge", data_type, gauge))


