library(cmdstanr)
library(posterior)
library(tidyverse)

mix_dens <- function(w, alphas, betas, weights) {
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}

mix_dens_sim <- function(w, alphas, betas, weights) {
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}

## LOGISITC ###
# filenames meaning:
# take3 = unordered, K = 6
# take2 = K=3, think it's unordered
# original = K=2, ordered
logistic_high_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "high_1_mix_betas_\\d{1}", full.names = TRUE)
post_param <- read_cmdstan_csv(logistic_high_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
logistic_high_1_data <- RcppSimdJson::fload("data/angular/logistic/high_1.json")
w1 <- logistic_high_1_data$w1
hist(w1, freq = FALSE, breaks = 60)
points(w1, mix_dens(w1, alphas, betas, weights), col = "blue")

logistic_high_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "high_1_mix_betas_take2_\\d{1}", full.names = TRUE)
post_param <- read_cmdstan_csv(logistic_high_1_files[2],variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)
alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
points(w1, mix_dens(w1, alphas, betas, weights), col = "orange")

logistic_high_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "high_1_mix_betas_take3_\\d{1}", full.names = TRUE)
post_param <- read_cmdstan_csv(logistic_high_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)
alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
points(w1, mix_dens(w1, alphas, betas, weights), col = "purple")

# Mid dependence
logistic_mid_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "mid_1_mix_betas_\\d{1}.csv", full.names = TRUE)
post_param <- read_cmdstan_csv(logistic_mid_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
logistic_mid_1_data <- RcppSimdJson::fload("data/angular/logistic/mid_1.json")
w1 <- logistic_mid_1_data$w1
hist(w1, freq = FALSE, breaks = 30)
points(w1, mix_dens(w1, alphas, betas, weights), col = "blue")

# low dependence
logistic_low_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "low_1_mix_betas_\\d{1}.csv", full.names = TRUE)
post_param <- read_cmdstan_csv(logistic_low_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
logistic_low_1_data <- RcppSimdJson::fload("data/angular/logistic/low_1.json")
w1 <- logistic_low_1_data$w1
hist(w1, freq = FALSE, breaks = 50)
points(w1, mix_dens(w1, alphas, betas, weights), col = "blue")

## GAUSS ##
# high dependence
gauss_high_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "high_1_mix_betas_\\d{1}.csv", full.names = TRUE)
post_param <- read_cmdstan_csv(gauss_high_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
gauss_high_1_data <- RcppSimdJson::fload("data/angular/gauss/high_1.json")
w1 <- gauss_high_1_data$w1
hist(w1, freq = FALSE, breaks = 50)
points(w1, mix_dens(w1, alphas, betas, weights), col = "blue")

gauss_high_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "high_1_mix_betas_3parts_unordered_\\d{1}.csv", full.names = TRUE)
post_param <- read_cmdstan_csv(gauss_high_1_files[1],variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
points(w1, mix_dens(w1, alphas, betas, weights), col = "orange")

# mid dependence
gauss_mid_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "mid_1_mix_betas_\\d{1}.csv", full.names = TRUE)
post_param <- read_cmdstan_csv(gauss_mid_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
gauss_mid_1_data <- RcppSimdJson::fload("data/angular/gauss/mid_1.json")
w1 <- gauss_mid_1_data$w1
hist(w1, freq = FALSE, breaks = 40)
points(w1, mix_dens(w1, alphas, betas, weights), col = "blue")

gauss_mid_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "mid_1_mix_betas_3parts_unordered_\\d{1}.csv", full.names = TRUE)
post_param <- read_cmdstan_csv(gauss_mid_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
points(w1, mix_dens(w1, alphas, betas, weights), col = "orange")

# low dependence
gauss_low_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "low_1_mix_betas_\\d{1}.csv", full.names = TRUE)
post_param <- read_cmdstan_csv(gauss_low_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
gauss_low_1_data <- RcppSimdJson::fload("data/angular/gauss/low_1.json")
w1 <- gauss_low_1_data$w1
hist(w1, freq = FALSE, breaks = 40)
points(w1, mix_dens(w1, alphas, betas, weights), col = "blue")

