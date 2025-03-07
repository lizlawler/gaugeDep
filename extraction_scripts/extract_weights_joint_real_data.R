args <- commandArgs(trailingOnly=TRUE)
likelihood <- args[1]
angle_dens <- args[2]

library(qs)
library(posterior)
library(loo)
library(dplyr)
library(tidyr)

data_type <- "redstone"

options(mc.cores = parallel::detectCores())

create_joint_loglik <- function(gauge, likelihood = "cens", angle_dens = "star") {
  trunc <- (likelihood == "trunc")
  star <- (angle_dens == "star")
  loglik_file <- sprintf("samplers/joint_loglik/real_data/%s_%s_%s_%s.qs", data_type, angle_dens, gauge, likelihood)
  
  # Check if joint loglikelihood file already exists
  if (file.exists(loglik_file)) {
    return(qread(loglik_file))
  }
  
  # If file doesn't already exist, create it
  temp_radial <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/real_data/pw_loglik/%s_%s_%s.qs", data_type, gauge, likelihood))
  
  temp_angular <- if (star) {
    qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/%s_%s.qs", data_type, gauge))
  } else {
    qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/pw_loglik/%s.qs", data_type))
  }
  
  if (trunc) {
    idx <- qread(sprintf("data/%s_expo.qs", data_type))$idx
    temp_angular <- temp_angular[, idx]
  }
  
  temp_joint <- temp_radial + temp_angular
  qsave(temp_joint, loglik_file)
  temp_joint
}

extract_lpd_pt <- function(gauge, likelihood, angle_dens) {
  temp <- create_joint_loglik(gauge = gauge, likelihood = likelihood, angle_dens = angle_dens)
  loo_temp <- loo(temp)
  return(loo_temp$pointwise[,"elpd_loo"])
}

create_lpd_list <- function(likelihood, angle_dens) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular") 
  lpd_list <- setNames(sapply(gauge_library,
                              function(x) extract_lpd_pt(gauge = x, likelihood = likelihood, angle_dens = angle_dens)),
                       gauge_library)
  return(lpd_list)
}

model_weights <- function(likelihood, angle_dens) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular") 
  temp <- create_lpd_list(likelihood = likelihood, 
                          angle_dens = angle_dens)
  stacking <- stacking_weights(temp)
  pseudobma_boot <- pseudobma_weights(temp)
  pseudobma_noboot <- pseudobma_weights(temp, BB = FALSE)
  wts <- tibble("stacking" = stacking, "pseudobma_boot" = pseudobma_boot, "pseudobma_noboot" = pseudobma_noboot) |>
    mutate(method = gauge_library,
           stacking = as.numeric(stacking),
           pseudobma_boot = as.numeric(pseudobma_boot),
           pseudobma_noboot = as.numeric(pseudobma_noboot))
  wts_file <- sprintf("fits_and_weights/wts_joint_model/%s_%s_%s.qs", data_type, likelihood, angle_dens)
  qsave(wts, wts_file)
  print(sprintf("Weights extracted and saved to disk for %s, with angular density: %s, and likelihood: %", angle_dens, likelihood))
}

model_weights(likelihood = likelihood, 
              angle_dens = angle_dens)
