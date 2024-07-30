library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
source("gauge_functions_wrt_x.R")

est_volume <- function(n = 500, pars = 0.5, gauge) {
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

mc_volume(N = 10000, pars = 0.1, "logistic")

# find true density of angles in L2-norm -------
polar_euc_tib <- function(x, y, gauge, dep) {
  euc_tib <- cbind(x,y) |> as_tibble()
  polar_euc_gw <- euc_tib |> mutate(r1 = x+y,
                                    r2 = sqrt(x^2 + y^2),
                                    w1 = x/r1,
                                    w2 = x/r2,
                                    gw1 = pmap_dbl(list(x = w1, y = 1-w1, dep = dep), function(x, y, dep) gauge(x, y, dep)),
                                    gw2 = pmap_dbl(list(x = w2, y = 1-w2, dep = dep), function(x, y, dep) gauge(x, y, dep)))
  return(polar_euc_gw)
}

## Sim data for varying dependence structures -----------
gauss <- function(N = 10000, dep = 0.5) {
  x <- rmvnorm(N, mean = c(0,0), sigma = matrix(c(1, dep, dep, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, gauss_gauge, dep))
}

logistic <- function(N = 10000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, logistic_gauge, dep))
}

inv_log <- function(N = 10000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log", mar1=c(1,1,1))
  y <- 1/x
  x <- y[,1]
  y <- y[,2]
  return(polar_euc_tib(x, y, inv_log_gauge, dep))
}

asym_log <- function(N = 10000, dep = 0.5, t1 = 0.5, t2 = 0.5) {
  x <- rbvevd(N, dep = dep, model = "alog", asy = c(t1, t2))
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, asym_log_gauge, dep))
}

dirichlet <- function(N = 10000, theta1, theta2 = 2) {
  x <- rbvevd(N, alpha = theta1, beta = theta2, model = 'ct')
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, dirichlet_gauge, dep))
}

test <- gauss(N = 10000, dep = 0.9)
w2_1 <- test$w2
w2_2 <- test$y / test$r2
w2_1norm <- w2_1 + w2_2

w1_1 <- test$w1
w1_2 <- test$y / test$r1
w1_2norm <- sqrt(w1_1^2 + w1_2^2)

w1_from_w2 <- w2_1 / w2_1norm
w2_from_w1 <- w1_1 / w1_2norm
hist(w2_from_w1)
hist(w2_1, add = TRUE, col = "blue")


test <- logistic(N = 10000, dep = 0.5)
r <- sqrt(test$x^2 + test$y^2)
w <- test$x / r
hist(acos(w), freq = FALSE, ylim = c(0,1.3))
# mc_volume(N = 10000, pars = 0.1, "logistic")
gw <- logistic_gauge(w, test$y/r, dep_par = 0.8)
test_dens <- (1/gw^2) / (2 * 0.8)
points(acos(w), test_dens, col = "blue")

r_test <- apply(test, 1, function(x) (sqrt(x[1]^2 + x[2]^2)))
w1_2norm <- test[,1] / r_test |> as_tibble()
ggplot(w1_2norm, aes(value)) + geom_histogram(binwidth = 0.05) + coord_polar()

plot(test)
test <- mc_volume(N=100, pars = 0.5, gauge = "gauss")
est_volume(n = 500, pars = c(0.5), "gauss")
n <- 500
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)
gx <- logistic_gauge(x1, x2, dep_par = 0.1)
sum(gx <= 1) / n

w <- seq(0, 1, length.out = 500)
gw <- logistic_gauge(w, (1-w), 0.9)
plot(w/gw, (1-w)/gw, type = "l", xlim = c(0,1))
which(((1-w)/gw) == 0)
((w)/gw)[500]
abline(v = 0.9, h = 0.9, col = "red")
est_volume(n = 500, pars = 0.5, "logistic")

test <- seq(0, 1.5, length.out=100)
sum(test <= 1)
