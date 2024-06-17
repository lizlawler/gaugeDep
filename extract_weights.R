args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
likelihood <- args[2]
threshold <- args[3]
level <- args[4]

library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(loo)

options(mc.cores = parallel::detectCores())

setwd("/data/accounts/lawler/research/gaugeDependence/")

create_model_fit <- function(sim_phase = "stacking", gauge, dep_type, likelihood, threshold, dep_level, dataset_num) {
  start_file_path <- paste0("stan/csv_fits/", sim_phase, "/", dep_type, "/", gauge, "/")
  csvfiles <- paste0(start_file_path,
                     list.files(path = start_file_path, 
                                pattern = paste0(dep_level, "_", dataset_num, "_", likelihood, "_", threshold, "_\\d{1}.csv")))
  fit <- read_cmdstan_csv(csvfiles, variables = "log_lik", format = "draws_matrix")$post_warmup_draws
  return(fit)
}

extract_lpd_pt <- function(sim_phase = "stacking", gauge, dep_type, likelihood, threshold, dep_level, dataset_num) {
  temp <- create_model_fit(gauge = gauge, dep_type = dep_type, likelihood = likelihood, threshold = threshold,
                           dep_level = dep_level, dataset_num = dataset_num)
  loo_temp <- loo(temp)
  return(loo_temp$pointwise[,"elpd_loo"])
}

create_lpd_list <- function(sim_phase = "stacking", dep_type, dep_level, likelihood, threshold, dataset_num) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular") 
  lpd_list <- setNames(sapply(gauge_library, 
                              function(x) extract_lpd_pt(gauge = x, dep_type = dep_type,
                                                         likelihood = likelihood, 
                                                         threshold = threshold,
                                                         dep_level = dep_level, 
                                                         dataset_num = dataset_num)), 
                       gauge_library)
  return(lpd_list)
}

model_weights <- function(sim_phase = "stacking", dep_type, dep_level, likelihood, threshold, dataset_num) {
  temp <- create_lpd_list(dep_type = dep_type, dep_level = dep_level, 
                          likelihood = likelihood, 
                          threshold = threshold,
                          dataset_num = dataset_num)
  stacking <- stacking_weights(temp)
  pseudobma_boot <- pseudobma_weights(temp)
  pseudobma_noboot <- pseudobma_weights(temp, BB = FALSE)
  print(paste0("Weights extracted for dataset: ", dataset_num))
  return(list("stacking" = stacking,
              "pseudobma_boot" = pseudobma_boot,
              "pseudobma_noboot" = pseudobma_noboot))
}

mod_wts <- lapply(1:100, function(x) model_weights(dep_type = dep_type, dep_level = level, likelihood = likelihood, 
                                                   threshold = threshold, dataset_num = x))
filepath <- paste0("stacking_weights/", dep_type, "_", level, "_", likelihood, "_", threshold, "_wts.RDS")
saveRDS(mod_wts, filepath)

print("Model weights have been successfully saved")