# Simulation of new pseudo-observations using the model structure
library(geometricMVE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(cubature)

csvfiles <- paste0("stan/csv_fits/calibrate/", "gauss", "/",
                   list.files(path = paste0("stan/csv_fits/calibrate/", "gauss", "/"), 
                              pattern = paste0("high", "_", 50, "_", "trunc", "_", "marg", "_\\d{1}.csv")))
fit <- as_cmdstan_fit(csvfiles)
med_pars <- fit |> as_draws_df() |> select(alpha, dep) |> apply(MARGIN = 2, FUN = median) |> as.numeric()
data <- fload("data/gauss/high_50.json")
r0_w <- data$r0_w
X1 <- data$R * data$W
X2 <- data$R * (1-data$W)
R <- data$R
W <- data$W

gauss_gauge <- function(x, y, rho = 0.5) {
  top <- x + y - 2 * rho * sqrt(x * y)
  return(top/(1-rho^2))
}

gw_fit <- gauss_gauge(W, 1-W, med_pars[2])
# gw_true <- gauss_gauge(W, 1-W, 0.9)
# plot(W/gw_fit, (1-W)/gw_fit, col = 4)
# points(W/gw_true, (1-W)/gw_true, col = 3)

rw_gw_fit <- cbind(R, W, r0_w) |> as_tibble() 
# ctau <- quantile(rw_gw_fit$gw_r, 0.99)
# rw_gw_fit <- rw_gw_fit |> mutate(r0w_tau = ctau / gw_fit)
rw_gw_fit_over1 <- rw_gw_fit |> filter(R > r0_w)

sim.2d_liz <- function (w, r0w, k = 1, nsim, par, gfun) 
{
  if (k != 1) {
    iw <- iweights.2d(k = k, r0w = r0w, w = w, gfun = gfun, 
                      par = par)
    star.ind <- sample(1:length(w), size = nsim, replace = T, 
                       prob = iw)
    wstar <- w[star.ind]
    r0w_star <- c(k * r0w)[star.ind]
  }
  else {
    star.ind <- sample(1:length(w), size = nsim, replace = T)
    wstar <- w[star.ind]
    r0w_star <- r0w[star.ind]
  }
  rate0 <- apply(cbind(w, 1 - w), 1, gfun, par = par[2:length(par)])
  rate <- rate0[star.ind]
  rstar <- qgamma(runif(nsim) * pgamma(r0w_star, shape = par[1], rate = rate, lower.tail = F) + pgamma(r0w_star, shape = par[1], rate = rate, lower.tail = T),
                  shape = par[1], rate = rate)
  xstar <- cbind(rstar * wstar, rstar * (1 - wstar))
  return(xstar)
}

# k = 1
new_x_liz <- sim.2d_liz(w=rw_gw_fit_over1$W, r0w=rw_gw_fit_over1$r0_w, k=1, 10000, par = med_pars, gfun = gauge_gaussian) |> as_tibble() |> rename(X1 = V1, X2=V2)
new_x <- sim.2d(w=rw_gw_fit_over1$W, r0w=rw_gw_fit_over1$r0_w, k=1, 10000, par = med_pars, gfun = gauge_gaussian) |> as_tibble() |> rename(X1 = V1, X2=V2)

plot(X1, X2,pch=20, xlim=c(0,14), ylim = c(0,14))
points(new_x_liz,pch=20,col=4)
points(new_x,pch=20,col=3)

lower <- 8
upper <- 10
prob_new_x_liz <- (new_x_liz |> filter(X1 >= lower & X1 <= upper & X2 >= lower & X2 <= upper) |> nrow())/nrow(new_x_liz) * nrow(rw_gw_fit_over1)/nrow(rw_gw_fit)
prob_new_x <- (new_x |> filter(X1 >= lower & X1 <= upper & X2 >= lower & X2 <= upper) |> nrow())/nrow(new_x) * nrow(rw_gw_fit_over1)/nrow(rw_gw_fit)

# convert [8,10] x [8,10] box in expo to gaussian coordinates
lower_limit <- qnorm(pexp(lower))
upper_limit <- qnorm(pexp(upper))
pmvnorm(lower = rep(lower_limit, 2), upper = rep(upper_limit, 2), corr = matrix(c(1, 0.9, 0.9, 1), nrow = 2))[1]


