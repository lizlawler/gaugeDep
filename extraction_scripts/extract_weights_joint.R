args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
likelihood <- args[3]

library(qs)
library(posterior)
library(loo)
library(dplyr)
library(tidyr)

options(mc.cores = parallel::detectCores())
angle_dens <- "vol"

create_joint_loglik <- function(dep_type, dep_level, gauge, data_num, likelihood = "cens", angle_dens = "vol") {
  trunc <- (likelihood == "trunc")
  vol <- (angle_dens == "vol")
  loglik_file <- sprintf("samplers/joint_loglik/%s/%s/%s_%s_%s_%s.qs",
                         dep_type, angle_dens, gauge, likelihood, dep_level, data_num)
  
  # Check if joint loglikelihood file already exists
  if (file.exists(loglik_file)) {
    return(qread(loglik_file))
  }
  
  # If file doesn't already exist, create it
  temp_radial <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/%s/pw_loglik/%s_%s_%s_%s.qs",
                               dep_type, gauge, likelihood, dep_level, data_num))
  
  temp_angular <- if (vol) {
    qread(sprintf("samplers/rcpp/angular_vol_mcmc_fits/%s/pw_loglik/%s_%s_%s.qs", 
                  dep_type, gauge, dep_level, data_num))
  } else {
    qread(sprintf("samplers/nimble/sb_mcmc_fits/%s/pw_loglik/%s_%s.qs", 
                  dep_type, dep_level, data_num))
  }
  
  if (trunc) {
    idx <- RcppSimdJson::fload(sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num))$idx
    temp_angular <- temp_angular[, idx]
  }
  
  temp_joint <- temp_radial + temp_angular
  qsave(temp_joint, loglik_file)
  temp_joint
}

extract_lpd_pt <- function(dep_type, dep_level, gauge, data_num, likelihood, angle_dens) {
  temp <- create_joint_loglik(dep_type = dep_type, dep_level = dep_level, gauge = gauge, 
                              data_num = data_num, likelihood = likelihood, angle_dens = angle_dens)
  loo_temp <- loo(temp)
  return(loo_temp$pointwise[,"elpd_loo"])
}

create_lpd_list <- function(dep_type, dep_level, data_num, likelihood, angle_dens) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular") 
  lpd_list <- setNames(sapply(gauge_library,
                              function(x) extract_lpd_pt(dep_type = dep_type, dep_level = dep_level, gauge = x,
                                                         data_num = data_num, likelihood = likelihood, angle_dens = angle_dens)),
                       gauge_library)
  return(lpd_list)
}

model_weights <- function(dep_type, dep_level, data_num, likelihood, angle_dens) {
  temp <- create_lpd_list(dep_type = dep_type, dep_level = dep_level, 
                          data_num = data_num,
                          likelihood = likelihood, 
                          angle_dens = angle_dens)
  stacking <- stacking_weights(temp)
  pseudobma_boot <- pseudobma_weights(temp)
  pseudobma_noboot <- pseudobma_weights(temp, BB = FALSE)
  print(paste0("Weights extracted for dataset: ", data_num))
  return(list("stacking" = stacking,
              "pseudobma_boot" = pseudobma_boot,
              "pseudobma_noboot" = pseudobma_noboot))
}

make_wts_df <- function(dep_type, dep_level, likelihood, angle_dens) {
  wts_file <- sprintf("fits_and_weights/wts_joint_model/%s_%s_%s_%s.qs",
                         dep_type, likelihood, angle_dens, dep_level)
  
  # Check if joint loglikelihood file already exists
  if (file.exists(wts_file)) {
    return(qread(wts_file))
  }
  
  wts <- lapply(1:100, function(x) model_weights(dep_type = dep_type, dep_level = dep_level, data_num = x,
                                                 likelihood = likelihood, angle_dens = angle_dens))
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  temp <- wts |>
    bind_rows() |> 
    mutate(method = rep(gauge_library, 100)) |>
    mutate(stacking = as.numeric(stacking),
           pseudobma_boot = as.numeric(pseudobma_boot),
           pseudobma_noboot = as.numeric(pseudobma_noboot)) |>
    mutate(dataset = rep(1:100, times = rep(6, 100)))

  qsave(temp, wts_file)
  temp
}

mod_wts <- make_wts_df(dep_type = dep_type, dep_level = dep_level, likelihood = likelihood, angle_dens = angle_dens)
print(sprintf("Model weights for %s, %s, %s, %s have been created and saved to disk", dep_type, dep_level, likelihood, angle_dens))
