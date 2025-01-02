args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
gauge <- args[3]
likelihood <- args[4]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

print(paste0("dep_type = ", dep_type))
print(paste0("dep_level = ", dep_level))
print(paste0("gauge = ", gauge))
print(paste0("likelihood = ", likelihood))

if(gauge != "dirichlet") {
  starting_vals <- c(rgamma(1, 4, 2), runif(1))
} else {
  starting_vals <- c(rgamma(1, 4, 2), abs(rt(1, 4,ncp = 0))*4, abs(rt(1, 4,ncp = 0))*2)
}

data_basename <- paste0("data/", dep_type, "/", dep_level, "_")

for(data_num in 1:100) {
  datafile <- paste0(data_basename, data_num, ".json")
  data <- RcppSimdJson::fload(datafile)
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
  qsave(x = results, file = paste0("samplers/rcpp/radial_mcmc_fits/", dep_type, "/", 
                                   gauge, "_", likelihood, "_", dep_level, "_", data_num, ".qs"))
  print(paste0("Successfully saved MCMC fit for dataset number: ", data_num))
}

