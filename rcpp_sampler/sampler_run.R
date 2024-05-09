args <- commandArgs(trailingOnly=TRUE)
gauge <- args[1]
likelihood <- args[2]
dataset_num <- args[3]

library(RcppSimdJson)
library(doParallel)
library(tidyverse)
Rcpp::sourceCpp("rcpp_sampler/gauge_mcmc.cpp")

file_path <- paste0("data/", gauge, "/", "mid_", dataset_num, ".json")
idx <- fload(file_path)$idx
W <- fload(file_path)$W
R <- fload(file_path)$R
r0_w <- fload(file_path)$r0_w

## logistic step_sizes
log_mid_steps <- c(0.2, 0.1)
# log_high_steps <- c(0.1, 0.01)
# log_low_steps <- c(0.2, 0.01)

## gauss step sizes
# gauss_low_steps <- c(0.3, 0.15)
gauss_mid_steps <- c(0.3, 0.15)
# gauss_high_steps <- c(0.15, 0.05)

gauss_mid_cens_steps <- c(0.1, 0.01)
log_mid_cens_steps <- c(0.1, 0.005)

# Obtain the current time with milliseconds
current_time <- format(Sys.time(), "%H%M%OS3")
# Convert time to a numeric value for seeding
current_seed <- as.numeric(gsub("[:.]", "", current_time))

set.seed(seed = current_seed, kind = "L'Ecuyer-CMRG")
registerDoParallel(5)
if(gauge == "gauss") {
  if(likelihood == "trunc") {
    results_mcmc <- foreach(c = 1:4) %dopar% {
      mcmc_mh(n_iter = 10000, W = W[idx], R = R[idx], r0_w = r0_w[idx], step_size = gauss_mid_steps, trunc = TRUE, gauss = TRUE)
    }
  } else {
    results_mcmc <- foreach(c = 1:4) %dopar% {
      mcmc_mh(n_iter = 5000, W = W, R = R, r0_w = r0_w, step_size = gauss_mid_cens_steps, trunc = FALSE, gauss = TRUE)
    }
  }
} else {
  if(likelihood == "trunc") {
    results_mcmc <- foreach(c = 1:4) %dopar% {
      mcmc_mh(n_iter = 10000, W = W[idx], R = R[idx], r0_w = r0_w[idx], step_size = log_mid_steps, trunc = TRUE, gauss = FALSE)
    }
  } else {
    results_mcmc <- foreach(c = 1:4) %dopar% {
      mcmc_mh(n_iter = 5000, W = W, R = R, r0_w = r0_w, step_size = log_mid_cens_steps, trunc = FALSE, gauss = FALSE)
    }
  }
}
stopImplicitCluster()

combined_chains <- map_df(results_mcmc, ~ as_tibble(.x, .name_repair = "unique"), .id = "chain")
names(combined_chains)[-1] <- c("iter", "alpha", "theta")
combined_chains <- combined_chains |> mutate(chain = as.factor(chain),
                                             iter = as.integer(iter))

# combined_chains |> ggplot(aes(x = iter, y = alpha, color = chain)) + geom_line()

write_csv(combined_chains, paste0("rcpp_sampler/csv_fits/", gauge, "/", likelihood, "_mid_", dataset_num, ".csv"))
