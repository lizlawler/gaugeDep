library(nimble)
library(tidyr)
library(dplyr)
library(qs)

data_type <- "redstone"

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

data <- qread(sprintf("data/%s_expo.qs", data_type))
params <- qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/%s.qs", data_type))
w <- data$W

results <- angular_loglik(angles = w, 
                          posterior_params = params)

qsave(x = results, file = sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/pw_loglik/%s.qs", data_type))
print(sprintf("Successfully saved posterior pointwise loglikelihood for angular mixture density of %s", data_type))
