args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
gauge <- args[3]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

print(paste0("dep_type = ", dep_type))
print(paste0("dep_level = ", dep_level))
print(paste0("gauge = ", gauge))

if(gauge != "dirichlet") {
  starting_vals <- runif(1)
} else {
  starting_vals <- c(abs(rt(1, 4,ncp = 0))*4, abs(rt(1, 4,ncp = 0))*2)
}

for(data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json",
                      dep_type, dep_level, data_num)
  data <- RcppSimdJson::fload(datafile)
  w <- data$W
  
  results <- angular_mcmc(angles = w, dim = 2, 
                          starting_theta = starting_vals, 
                          gauge_type = gauge, 
                          n_updates = 15000, 
                          update_freq = 250, 
                          n_burnin = 5000, 
                          n_thin = 5, 
                          adapt_cov = TRUE)
  qsave(x = results, file = sprintf("samplers/rcpp/angular_vol_mcmc_fits/%s/%s_%s_%s.qs",
                                    dep_type, gauge, dep_level, data_num))
  print(paste0("Successfully saved MCMC fit for dataset number: ", data_num))
}

