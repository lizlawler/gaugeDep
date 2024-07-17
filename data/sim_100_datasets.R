library(evd)
library(mvtnorm)
library(tidyverse)

polar_euc_tib <- function(x, y, gauge, dep) {
  euc_tib <- cbind(x,y) |> as_tibble()
  polar_euc_gw <- euc_tib |> mutate(r = x+y,
                                 w = x/r,
                                 gw = pmap_dbl(list(x = w, y = 1-w, dep = dep), function(x, y, dep) gauge(x, y, dep)),
                                 gw_r = pmap_dbl(list(x = gw, y = r), function(x, y) x*y))
  return(polar_euc_gw)
}

# for Husler-Reiss data only (doesn't use gauge function threshold)
polar_euc_tib_hr <- function(x, y) {
  euc_tib <- cbind(x,y) |> as_tibble()
  polar_euc <- euc_tib |> mutate(r = x+y,
                                    w = x/r)
  return(polar_euc)
}

## 0.25, 2, 6

## Copula functions for varying dependence structures -----------
gauss_gauge <- function(x, y, dep = 0.5) {
  top <- x + y - 2 * dep * sqrt(x * y)
  return(top/(1-dep^2))
}

gauss <- function(N = 10000, dep = 0.5) {
  x <- rmvnorm(N, mean = c(0,0), sigma = matrix(c(1, dep, dep, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, gauss_gauge, dep))
}

logistic_gauge <- function(x, y, dep = 0.5) {
  dep_inv <- 1/dep
  return(dep_inv * pmax(x, y) + (1-dep_inv) * pmin(x,y))
}

logistic <- function(N = 10000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, logistic_gauge, dep))
}

inv_log_gauge <- function(x, y, dep = 0.5) ((x^(1/dep) + y^(1/dep))^dep)

inv_log <- function(N = 10000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log", mar1=c(1,1,1))
  y <- 1/x
  x <- y[,1]
  y <- y[,2]
  return(polar_euc_tib(x, y, inv_log_gauge, dep))
}

asym_log_gauge <- function(x, y, dep = 0.5) {
  dep_inv <- 1/dep
  return(pmin((x + y), (dep_inv * pmax(x, y) + (1-dep_inv)*pmin(x,y))))
}

asym_log <- function(N = 10000, dep = 0.5, t1 = 0.5, t2 = 0.5) {
  x <- rbvevd(N, dep = dep, model = "alog", asy = c(t1, t2))
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, asym_log_gauge, dep))
}

dirichlet_gauge <- function(x, y, theta1, theta2 = 2) {
  return((1 + theta1 + theta2) * pmax(x, y) - (theta1 * x + theta2 * y))
}

dirichlet <- function(N = 10000, theta1, theta2 = 2) {
  x <- rbvevd(N, alpha = theta1, beta = theta2, model = 'ct')
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, dirichlet_gauge, dep))
}

husler_reiss <- function(N = 10000, dep = 1) {
  x <- rbvevd(N, dep = dep, model = "hr")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib_hr(x, y))
}

## Creation of stan data list -------------
grab_top_n <- function(cloud_tib, n0 = 1, N = 10000) {
  tau <- (N-n0)/N
  q1 <- quantile(cloud_tib$x, tau)
  q2 <- quantile(cloud_tib$y, tau)
  q <- max(q1, q2)
  idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  eps <- 0.001
  while (length(idx) > n0) {
    q <- q + eps
    idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  }
  r0_w <- cloud_tib |> mutate(r0_w = ifelse(w > 0.5, q/w, q/(1-w))) |> select(r0_w)
  if("gw" %in% colnames(cloud_tib)) {
    ctau <- as.numeric(quantile(cloud_tib$gw_r, tau))
    r0_w_ctau <- ctau / cloud_tib$gw
    idx_ctau <- which(cloud_tib$r > r0_w_ctau)
    return(list(q = q,
                idx = idx,
                n0 = length(idx),
                N = N,
                R = cloud_tib$r,
                W = cloud_tib$w,
                r0_w = r0_w$r0_w,
                ctau = ctau,
                r0_w_ctau = r0_w_ctau,
                n0_ctau = length(idx_ctau),
                idx_ctau = idx_ctau))
  } else {
    return(list(q = q,
                idx = idx,
                n0 = length(idx),
                N = N,
                R = cloud_tib$r,
                W = cloud_tib$w,
                r0_w = r0_w$r0_w))
  }
}


data_list <- function(N = 10000, n0 = 1, dep, cop_func) {
  og_data <- cop_func(N, dep)
  return(grab_top_n(og_data, n0 = n0, N = N))
}

gen_data_file <- function(N = 10000, n0 = 1, dep, cop_func, dep_level, iter) {
  temp <- data_list(N, n0, dep, cop_func)
  cmdstanr::write_stan_json(temp, paste0("data/", deparse(substitute(cop_func)), "/", dep_level, "_", iter, ".json"))
  # cmdstanr::write_stan_json(temp, paste0("data/", dep_type, "/", dep_level, "_", iter, ".json"))
}

# test <- data_list(dep = 6, cop_func = husler_reiss)
dep_levels <- list(c(0.25, "low"), c(2, "mid"), c(6, "high"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(5000, 250, dep, husler_reiss, level, i)
  }
}


dep_levels <- list(c(0.1, "low"), c(0.5, "mid"), c(0.9, "high"), c(0.8, "wc"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(5000, 250, dep, gauss, level, i)
  }
}

data <- RcppSimdJson::fload("data/gauss/wc_10.json")
x1 <- data$R * data$W 
x2 <- data$R * (1 - data$W)
idx <- data$idx
idx_ctau <- data$idx_ctau
plot(x1, x2, pch = 20)
points(x1[idx], x2[idx], pch = 20, col = "blue")
points(x1[idx_ctau], x2[idx_ctau], pch = 20, col = "red")

dep_levels <- list(c(0.1, "high"), c(0.5, "mid"), c(0.9, "low"), c(0.8, "wc_low"), c(0.4, "wc_mid"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(5000, 250, dep, logistic, level, i)
  }
}

data <- RcppSimdJson::fload("data/logistic/high_10.json")
x1 <- data$R * data$W 
x2 <- data$R * (1 - data$W)
idx <- data$idx
idx_ctau <- data$idx_ctau
plot(x1, x2, pch = 20)
points(x1[idx], x2[idx], pch = 20, col = "blue")
points(x1[idx_ctau], x2[idx_ctau], pch = 20, col = "red")

dep_levels <- list(c(3, "high"), c(1, "mid"), c(0.5, "low"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:100) {
    gen_data_file(10000, 1000, dep, dirichlet, level, i)
  }
}
