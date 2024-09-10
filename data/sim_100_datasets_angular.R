library(evd)
library(mvtnorm)
library(tidyverse)
library(cmdstanr)
source("gauge_functions_wrt_x.R")
# 
model <- cmdstan_model("stan/radial_angular/bivar_cens_marg_logistic_angular.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)
model$expose_functions(global = TRUE)
# MC volume estimation
est_volume <- function(n = 1000, pars = 0.5, gauge) {
  x1 <- runif(n, 0, 1)
  x2 <- runif(n, 0, 1)
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  gx <- gauge_fcn(x1, x2, dep_par = pars)
  return(sum(gx <= 1) / n)
}

mc_volume <- function(N = 10000, pars = 0.5, gauge) {
  rep_est <- replicate(N, est_volume(n = 500, pars = pars, gauge = gauge))
  return(mean(rep_est))
}

temp <- RcppSimdJson::fload("data/angular/logistic/high_1.json")
temp$x2 <- grid[,2]
temp$n_grid <- nrow(grid)
write_stan_json(temp, "data/angular/logistic/high_1.json")

x1_vec <- seq(0, 1, length.out=100)
grid <- expand.grid(x1_vec, x1_vec)
## Copula functions for varying dependence structures -----------
gauss <- function(N = 10000, dep = 0.5) {
  x <- rmvnorm(N, mean = c(0,0), sigma = matrix(c(1, dep, dep, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  r <- x + y
  w1 <- x/r
  return(cbind(r, w1) |> as_tibble())
}

logistic <- function(N = 10000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  r <- x + y
  w1 <- x/r
  return(cbind(r, w1) |> as_tibble())
}


## Creation of stan data list -------------
data_list <- function(N = 10000, dep, cop_func) {
  if(cop_func == "logistic") {
    L <- dep
  } else {
    L <- mc_volume(10000, pars = dep, gauge = cop_func)
  }
  cop_func <- get(cop_func)
  data_tib <- cop_func(N, dep)
  return(list(N = N,
              d = 2,
              w1 = data_tib$w1,
              r = data_tib$r,
              L = L))
}

gen_data_file <- function(N = 10000, dep, cop_func, dep_level, iter) {
  temp <- data_list(N, dep, cop_func)
  cmdstanr::write_stan_json(temp, paste0("data/angular/", cop_func, "/", dep_level, "_", iter, ".json"))
}

dep_levels <- list(c(0.1, "low"), c(0.5, "mid"), c(0.9, "high"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(5000, dep, "gauss", level, i)
  }
}

x1 <- seq(0, 1, length.out = 100)
grid <- expand.grid(x1, x1)
x1 <- grid[,1]
x2 <- grid[,2]
n_grid <- as.numeric(nrow(grid))
# read in previously made datasets to append grid values for estimating volume with the dependence parameter in Stan model
append_data <- function(level, data_num) {
  file <- paste0("data/gauss/", level, "_", data_num, ".json")
  temp <- RcppSimdJson::fload(file)
  temp$x1 <- x1
  temp$x2 <- x2
  temp$n_grid <- n_grid
  cmdstanr::write_stan_json(temp, file)
}

gauss_data_combos <- expand_grid(levels = c("low", "mid", "high"), 
                                 datasets = as.character(1:100))
apply(gauss_data_combos, 1, function(row) {
  append_data(row["levels"], row["datasets"])})

k_vals <- list(gauss = list(low = 5, mid = 5, high = 19),
               logistic = list(low = 5, mid = 9, high = 89))
# read in previously made datasets to append K values for estimating density with Bernstein polynomial
append_data <- function(dep_type, level, data_num) {
  file <- paste0("data/angular/", dep_type, "/", level, "_", data_num, ".json")
  temp <- RcppSimdJson::fload(file)
  temp$K <- k_vals[[dep_type]][[level]]
  cmdstanr::write_stan_json(temp, file)
}

types_data_combos <- expand_grid(dep_types = c("gauss", "logistic"),
                                 levels = c("low", "mid", "high"), 
                                 datasets = as.character(1:100))
apply(types_data_combos, 1, function(row) {
  append_data(row["dep_types"],
              row["levels"], 
              row["datasets"])})

dep_levels <- list(c(0.1, "high"), c(0.5, "mid"), c(0.9, "low"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(5000, dep, "logistic", level, i)
  }
}
