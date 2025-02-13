args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
gauge <- args[3]
likelihood <- args[4]

library(gaugeDependence)
library(dplyr)
library(tidyr)
library(qs)

for(data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json",
                      dep_type, dep_level, data_num)
  paramsfile <- sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs",
                        dep_type, gauge, likelihood, dep_level, data_num)
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
  qsave(x = results, file = sprintf("samplers/rcpp/radial_mcmc_fits/%s/pw_loglik/%s_%s_%s_%s.qs",
                                    dep_type, gauge, likelihood, dep_level, data_num))
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}