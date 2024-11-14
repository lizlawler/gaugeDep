library(tidyverse)
library(doParallel)
library(foreach)
library(gaugeDep)

data <- RcppSimdJson::fload("data/gauss/high_1.json")
W <- data$W
R <- data$R
r0w <- data$r0_w
grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))
sum_grid <- grid[,1] + grid[,2]
sqrt_grid <- sqrt(grid[,1] * grid[,2])
results <- radial_angular_mcmc(R, W, r0w, 
                    sum_term = sum_grid, sqrt_term = sqrt_grid, dim = 2,
                    starting_theta_rad = c(rgamma(1, 4, 2), runif(1)), starting_theta_ang = runif(1), 
                    rad_gauge_type = "gauss", ang_gauge_type = "gauss",
                    n_updates = 15000, update_freq = 250,
                    adapt_cov = TRUE)

plot(results$samples[,"alpha"][5000:15000], type = "l")

RNGkind("L'Ecuyer-CMRG")
registerDoParallel(8)

gauss_samples <-
  foreach(levels = c("low", "mid", "high")) %:%
    foreach(dataset = 1:100) %:%
      foreach(chain = 1:3) %dopar% {
        data <- RcppSimdJson::fload(paste0("data/gauss/", levels, "_", dataset, ".json"))
        W <- data$W
        R <- data$R
        r0w <- data$r0_w
        data_list <- list(R = R, W = W)
        grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))
        sum_grid <- grid[,1] + grid[,2]
        sqrt_grid <- sqrt(grid[,1] * grid[,2])
        radial_angular_mcmc(R, W, r0w, 
                            sum_term = sum_grid, sqrt_term = sqrt_grid, dim = 2,
                            starting_theta_rad = c(rgamma(1, 4, 2), runif(1)), starting_theta_ang = runif(1), 
                            rad_gauge_type = "gauss", ang_gauge_type = "gauss",
                            n_updates = 15000, update_freq = 250,
                            adapt_cov = TRUE)
        
      }
stopImplicitCluster()

walk2(gauss_samples, c("low", "mid", "high"), function(level_samples, level_name) {
  walk2(level_samples, seq_along(level_samples), function(dataset_samples, dataset_num) {
    named_list <- set_names(dataset_samples, c("chain1", "chain2", "chain3"))
    saveRDS(named_list, file = paste0("mcmc_samples/gauss/", level_name, "_", dataset_num, ".rds"))
  })
})

RNGkind("L'Ecuyer-CMRG")
registerDoParallel(8)

logistic_samples <-
  foreach(levels = c("low", "mid", "high")) %:%
    foreach(dataset = 1:100) %:%
      foreach(chain = 1:3) %dopar% {
        data <- RcppSimdJson::fload(paste0("data/logistic/", levels, "_", dataset, ".json"))
        W <- data$W
        R <- data$R
        r0w <- data$r0_w
        data_list <- list(R = R, W = W)
        grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))
        sum_grid <- grid[,1] + grid[,2]
        sqrt_grid <- sqrt(grid[,1] * grid[,2])
        radial_angular_mcmc(R, W, r0w, 
                            sum_term = sum_grid, sqrt_term = sqrt_grid, dim = 2,
                            starting_theta_rad = c(rgamma(1, 4, 2), runif(1)), starting_theta_ang = runif(1), 
                            rad_gauge_type = "logistic", ang_gauge_type = "gauss",
                            n_updates = 15000, update_freq = 250,
                            adapt_cov = TRUE)
        
      }
stopImplicitCluster()

walk2(logistic_samples, c("low", "mid", "high"), function(level_samples, level_name) {
  walk2(level_samples, seq_along(level_samples), function(dataset_samples, dataset_num) {
    named_list <- set_names(dataset_samples, c("chain1", "chain2", "chain3"))
    saveRDS(named_list, file = paste0("mcmc_samples/logistic/", level_name, "_", dataset_num, ".rds"))
  })
})


# ## reshape logistic rds files
# reshape_rds_file <- function(file_path) {
#   # Load the .rds file
#   mcmc_chains <- readRDS(file_path)
#
#   params_list <- lapply(mcmc_chains, function(chain) {
#     alpha <- chain$trace[,1]
#     dep_w <- chain$trace[,2]
#     dep_r <- chain$trace[,3]
#     return(list(alpha = alpha, dep_w = dep_w, dep_r = dep_r,
#                 sigma_m_trace = chain$sigma_m_trace,
#                 r_trace = chain$r_trace,
#                 acc_prob = chain$acc_prob))
#   })
#
#   # Save the modified data back into the file (or a new file)
#   saveRDS(params_list, file = file_path)
# }
# directory <- "mcmc_samples/logistic/"
# # Get all the .rds file paths in the directory
# rds_files <- list.files(directory, pattern = "\\.rds$", full.names = TRUE)
#
# # Use purrr::walk to apply the function to all files
# purrr::walk(rds_files, reshape_rds_file)
