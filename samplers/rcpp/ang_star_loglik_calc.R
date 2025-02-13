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

for(data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  paramsfile <- sprintf("samplers/rcpp/angular_vol_mcmc_fits/%s/%s_%s_%s.qs",
                        dep_type, gauge, dep_level, data_num)
  data <- RcppSimdJson::fload(datafile)
  params <- qread(paramsfile)$samples
  params <- params[,1:(ncol(params) - 1)]
  w <- data$W
  
  results <- angular_loglik(angles = w, dim = 2,
                            posterior_params = params, 
                            gauge_type = gauge)
  
  savename <- sprintf("samplers/rcpp/angular_vol_mcmc_fits/%s/pw_loglik/%s_%s_%s.qs",
                      dep_type, gauge, dep_level, data_num)
  qsave(x = results, file = savename)
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}

