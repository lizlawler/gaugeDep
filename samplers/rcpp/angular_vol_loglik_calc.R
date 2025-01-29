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

data_basename <- paste0("data/", dep_type, "/", dep_level, "_")
params_basename <- paste0("samplers/rcpp/angular_vol_mcmc_fits/", 
                          dep_type, "/",
                          gauge, "_", dep_level, "_")

for(data_num in 1:100) {
  datafile <- paste0(data_basename, data_num, ".json")
  paramsfile <- paste0(params_basename, data_num, ".qs")
  data <- RcppSimdJson::fload(datafile)
  params <- qread(paramsfile)$samples
  params <- params[,1:(ncol(params) - 1)]
  w <- data$W
  
  results <- angular_loglik(angles = w, dim = 2,
                          posterior_params = params, 
                          gauge_type = gauge)
  qsave(x = results, file = paste0("samplers/rcpp/angular_vol_mcmc_fits/", dep_type, "/pw_loglik/", 
                                   gauge, "_", dep_level, "_", data_num, ".qs"))
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}

