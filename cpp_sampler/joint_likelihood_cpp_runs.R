library(tidyverse)
library(doParallel)
library(foreach)
Rcpp::sourceCpp("radial_angular_mcmc.cpp")

RNGkind("L'Ecuyer-CMRG")
registerDoParallel(6)

samples <- 
  foreach(dataset = 1:100) %:%
    foreach(chain = 1:3) %dopar% {
      data <- RcppSimdJson::fload(paste0("data/gauss/high_", 1, ".json"))
      W <- data$W
      R <- data$R
      r0w <- data$r0_w
      data_list <- list(R = R, W = W)
      grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))
      sum_grid <- grid[,1] + grid[,2]
      sqrt_grid <- sqrt(grid[,1] * grid[,2])
      adaptive_mh_rad_ang_vol(z = data_list,
                              threshold = r0w,
                              sum_term = sum_grid,
                              sqrt_term = sqrt_grid,
                              dim = 2,
                              starting_theta = c(rgamma(1, 4, 2), runif(1), runif(1)),
                              n_updates = 25000, update_freq = 500,
                              adapt_cov = TRUE)
      
    }
stopImplicitCluster()

library(purrr)
walk2(samples, seq_along(samples), ~ {
  named_list <- set_names(.x, c("chain1", "chain2", "chain3"))
  saveRDS(named_list, file = paste0("mcmc_samples/gauss/high_", .y, ".rds"))
})
