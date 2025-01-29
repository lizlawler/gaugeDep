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

data_basename <- paste0("data/", dep_type, "/", dep_level, "_")
params_basename <- paste0("samplers/rcpp/radial_mcmc_fits/", 
                          dep_type, "/",
                          gauge, "_", likelihood, "_", dep_level, "_")

for(data_num in 1:100) {
  datafile <- paste0(data_basename, data_num, ".json")
  paramsfile <- paste0(params_basename, data_num, ".qs")
  data <- RcppSimdJson::fload(datafile)
  params <- qread(paramsfile)$samples
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
  qsave(x = results, file = paste0("samplers/rcpp/radial_mcmc_fits/", dep_type, "/pw_loglik/",
                                   gauge, "_", likelihood, "_", dep_level, "_", data_num, ".qs"))
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}