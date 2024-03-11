library(evd)
library(mvtnorm)
library(tidyverse)


## Copula functions for varying dependence structures -----------
gauss <- function(n = 1000, r = 0.5) {
  x <- rmvnorm(n, mean = c(0,0), sigma = matrix(c(1, r, r, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y))
}

logistic <- function(n = 1000, r = 0.5) {
  x <- rbvevd(n, dep = r, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y))
}

inv_log <- function(n = 1000, r = 0.5) {
  x <- rbvevd(n, dep = r, model = "log", mar1=c(1,1,1))
  y <- 1/x
  x <- y[,1]
  y <- y[,2]
  return(cbind(x,y))
}

asym_log <- function(n = 1000, r = 0.5, t1 = 0.5, t2 = 0.5) {
  x <- rbvevd(n, dep = r, model = "alog", asy = c(t1, t2))
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y))
}

## Creation of stan data list -------------
data_list <- function(n, dep, cop_func) {
  og_data <- cop_func(n, dep)
  trunc_data <- og_data %>% as_tibble() %>%
    mutate(q1 = quantile(og_data[,1], 0.95), 
           q2 = quantile(og_data[,2], 0.95),
           high = case_when(x >= q1 | y >= q2 ~ 1,
                            .default = 0),
           high = as.factor(high),
           r = x + y,
           w = x/r,
           r0_w = ifelse(w > 0.5, q1/w, q2/(1-w))) %>% 
    filter(high == 1) %>% select(r, w, r0_w)
  stan_data <- list(R = trunc_data$r, 
                    W = trunc_data$w,
                    r0_w = trunc_data$r0_w,
                    N = nrow(trunc_data))
  return(stan_data)
}

gen_data_file <- function(n, dep, cop_func, dep_level, iter) {
  temp <- data_list(n, dep, cop_func)
  cmdstanr::write_stan_json(temp, paste0("data/", deparse(substitute(cop_func)), "/", dep_level, "_", iter, ".json"))
}

dep_levels <- list(c(0.1, "low"), c(0.5, "mid"), c(0.9, "high"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(10000, dep, gauss, level, i)
  }
}

dep_levels <- list(c(0.1, "high"), c(0.5, "mid"), c(0.9, "low"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(10000, dep, logistic, level, i)
  }
}

