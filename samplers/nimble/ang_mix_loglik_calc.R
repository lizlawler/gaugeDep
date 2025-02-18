args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]

library(nimble)
library(tidyr)
library(dplyr)
library(qs)

mix_lpdf <- function(w, wts, alpha, beta) {
  n <- length(wts)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + wts[i] * dbeta(w, alpha[i], beta[i])
  }
  return(dens)
}

angular_loglik <- function(angles, posterior_params) {
  probs <- posterior_params[,grepl(pattern = "probs\\[", colnames(posterior_params))]
  alpha <- posterior_params[,grepl(pattern = "alphastar\\[", colnames(posterior_params))]
  beta <- posterior_params[,grepl(pattern = "betastar\\[", colnames(posterior_params))]
  
  n_iter <- nrow(probs)
  n_obs <- length(angles)
  
  pw_loglik <- matrix(NA, nrow = n_iter, ncol = n_obs)
  for(i in 1:n_iter) {
    pw_loglik[i, ] <- mix_lpdf(angles, probs[i,], alpha[i,], beta[i,])
  }
  return(pw_loglik) 
}

for(data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
  paramsfile <- sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
                        dep_type, dep_level, data_num)
  
  data <- RcppSimdJson::fload(datafile)
  params <- qread(paramsfile)
  w <- data$W
  
  results <- angular_loglik(angles = w, 
                            posterior_params = params)
  savename <- sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/pw_loglik/%s_%s.qs",
                      dep_type, dep_level, data_num)
  qsave(x = results, file = savename)
  print(paste0("Successfully saved posterior pointwise loglikelihood for dataset number: ", data_num))
}
