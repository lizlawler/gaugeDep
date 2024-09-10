library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
source("gauge_functions_wrt_x.R")

est_volume <- function(n = 200, pars = 0.5, gauge) {
  # x1 <- runif(n, 0, 1)
  # x2 <- runif(n, 0, 1)
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  gx <- gauge_fcn(grid[,1], grid[,2], dep_par = pars)
  return(sum(gx <= 1) / n^2)
}

n <- 100
temp <- seq(0, 1, length.out = n)
grid <- expand.grid(temp, temp)
gx <- logistic_gauge(grid[,1], grid[,2], dep_par = 0.9)

inv_logit_test <- function(x) {
  return(1 / (1 + exp(-x)))
}

temp_vec <- approx_indicator(gx, 15)
approx_indicator <- function(x, k) {
  temp_arg <- k * (1.0 - x)
  return(inv_logit_test(temp_arg))
}
x_vals <- seq(-3, 3, length.out = 5000)
plot(x_vals, approx_indicator(x_vals, 10), type = "l", col = "red")
lines(x_vals, approx_indicator(x_vals, 15), col = "blue")
lines(x_vals, approx_indicator(x_vals, 50), col = "orange")

indic_vec <- rep(NA, n)
for (i in 1:n^2) {
    indic_vec[i] = inv_logit_test(-500 * (gx[i] - 1.0));
}
mean(indic_vec)

inv_logit_test(-10*(gx[1] - 1.0))
# real est_vol(int n_vol, vector sum_term, vector sqrt_term, real pars) {
#   vector[n_vol] gx = (sum_term - 2 * pars * sqrt_term) / (1 - pars^2);
#   vector[n_vol] logistic_indicator_vol;
#   real k = 10; // You can adjust the steepness parameter 'k' as needed.
#   for (i in 1:n_vol) {
#     logistic_indicator_vol[i] = inv_logit(-k * (gx[i] - 1.0));
#   }
#   return mean(logistic_indicator_vol);
# }


mc_volume <- function(N = 10000, pars = 0.5, gauge) {
  rep_est <- replicate(N, est_volume(n = 100, pars = pars, gauge = gauge))
  return(mean(rep_est))
}

## Sim data for varying dependence structures -----------
gauss <- function(N = 10000, dep = 0.5) {
  x <- rmvnorm(N, mean = c(0,0), sigma = matrix(c(1, dep, dep, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y) |> as_tibble())
}

logistic <- function(N = 10000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y) |> as_tibble())
}

inv_log <- function(N = 10000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log", mar1=c(1,1,1))
  y <- 1/x
  x <- y[,1]
  y <- y[,2]
  return(cbind(x, y) |> as_tibble())
}

asym_log <- function(N = 10000, dep = 0.5, t1 = 0.5, t2 = 0.5) {
  x <- rbvevd(N, dep = dep, model = "alog", asy = c(t1, t2))
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y) |> as_tibble())
}

dirichlet <- function(N = 10000, dep_par) {
  theta1 <- dep_par[1]
  theta2 <- dep_par[2]
  x <- rbvevd(N, alpha = theta1, beta = theta2, model = 'ct')
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y) |> as_tibble())
}

# density functions ---------
# spherical angles
dens_sphere <- function(w1, w2, gauge, par_val, norm_type) {
  if(gauge == "logistic") {
    mc_vol <- par_val
  } else {
    mc_vol <- mc_volume(N = 10000, pars = par_val, gauge)
  }
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  gw <- gauge_fcn(w1, w2, par_val)
  if(norm_type == "L2") {
    return(1 / (gw^2 * 2 * mc_vol))
  } else {
    return((w1_1^2 + w1_2^2)/ (gw^2 * 2 * mc_vol))
  }
}

# L1, pseudo angles
dens_l1_norm <- function(w1, w2, gauge, par_val) {
  if(gauge == "logistic") {
    mc_vol <- par_val
  } else {
    mc_vol <- mc_volume(N = 10000, pars = par_val, gauge)
  }
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  gw <- gauge_fcn(w1, w2, par_val)
  return(1 / (gw^2 * 2 * mc_vol))
}

# Logistic cases ------
# high - spherical angle density
dep_val <- 0.1
logistic_high <- logistic(N = 10000, dep = dep_val)
r2 <- sqrt(logistic_high$x^2 + logistic_high$y^2)
w2_1 <- logistic_high$x / r2
w2_2 <- logistic_high$y / r2
r1 <- logistic_high$x + logistic_high$y
w1_1 <- logistic_high$x / r1
w1_2 <- logistic_high$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "logistic", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "logistic", dep_val, "L1"), col = "green")

# density in terms of L1 norm and L1 angles
hist(w1_1, freq = FALSE, ylim = c(0,15))
points(w1_1, dens_l1_norm(w1_1, w1_2, "logistic", dep_val))

# logistic - mid ------
dep_val <- 0.5
logistic_mid <- logistic(N = 10000, dep = dep_val)
r2 <- sqrt(logistic_mid$x^2 + logistic_mid$y^2)
w2_1 <- logistic_mid$x / r2
w2_2 <- logistic_mid$y / r2
r1 <- logistic_mid$x + logistic_mid$y
w1_1 <- logistic_mid$x / r1
w1_2 <- logistic_mid$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "logistic", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "logistic", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "logistic", dep_val))

# logistic - low ------
dep_val <- 0.9
logistic_low <- logistic(N = 10000, dep = dep_val)
r2 <- sqrt(logistic_low$x^2 + logistic_low$y^2)
w2_1 <- logistic_low$x / r2
w2_2 <- logistic_low$y / r2
r1 <- logistic_low$x + logistic_low$y
w1_1 <- logistic_low$x / r1
w1_2 <- logistic_low$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "logistic", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "logistic", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "logistic", dep_val))

# gauss cases ------
# high - spherical angle density (using different norms but still spherical angles)
dep_val <- 0.9
gauss_high <- gauss(N = 10000, dep = dep_val)
r2 <- sqrt(gauss_high$x^2 + gauss_high$y^2)
w2_1 <- gauss_high$x / r2
w2_2 <- gauss_high$y / r2
r1 <- gauss_high$x + gauss_high$y
w1_1 <- gauss_high$x / r1
w1_2 <- gauss_high$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "gauss", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "gauss", dep_val, "L1"), col = "green")

# density in terms of L1 norm and L1 angles
hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "gauss", dep_val))

# gauss - mid ------
dep_val <- 0.5
gauss_mid <- gauss(N = 10000, dep = dep_val)
r2 <- sqrt(gauss_mid$x^2 + gauss_mid$y^2)
w2_1 <- gauss_mid$x / r2
w2_2 <- gauss_mid$y / r2
r1 <- gauss_mid$x + gauss_mid$y
w1_1 <- gauss_mid$x / r1
w1_2 <- gauss_mid$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "gauss", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "gauss", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "gauss", dep_val))

# gauss - low ------
dep_val <- 0.1
gauss_low <- gauss(N = 10000, dep = dep_val)
r2 <- sqrt(gauss_low$x^2 + gauss_low$y^2)
w2_1 <- gauss_low$x / r2
w2_2 <- gauss_low$y / r2
r1 <- gauss_low$x + gauss_low$y
w1_1 <- gauss_low$x / r1
w1_2 <- gauss_low$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "gauss", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "gauss", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "gauss", dep_val))

# inv_log cases ------
# high - spherical angle density (using different norms but still spherical angles)
dep_val <- 0.1
inv_log_high <- inv_log(N = 10000, dep = dep_val)
r2 <- sqrt(inv_log_high$x^2 + inv_log_high$y^2)
w2_1 <- inv_log_high$x / r2
w2_2 <- inv_log_high$y / r2
r1 <- inv_log_high$x + inv_log_high$y
w1_1 <- inv_log_high$x / r1
w1_2 <- inv_log_high$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "inv_log", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "inv_log", dep_val, "L1"), col = "green")

# density in terms of L1 norm and L1 angles
hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "inv_log", dep_val))

# inv_log - mid ------
dep_val <- 0.5
inv_log_mid <- inv_log(N = 10000, dep = dep_val)
r2 <- sqrt(inv_log_mid$x^2 + inv_log_mid$y^2)
w2_1 <- inv_log_mid$x / r2
w2_2 <- inv_log_mid$y / r2
r1 <- inv_log_mid$x + inv_log_mid$y
w1_1 <- inv_log_mid$x / r1
w1_2 <- inv_log_mid$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "inv_log", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "inv_log", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "inv_log", dep_val))

# inv_log - low ------
dep_val <- 0.9
inv_log_low <- inv_log(N = 10000, dep = dep_val)
r2 <- sqrt(inv_log_low$x^2 + inv_log_low$y^2)
w2_1 <- inv_log_low$x / r2
w2_2 <- inv_log_low$y / r2
r1 <- inv_log_low$x + inv_log_low$y
w1_1 <- inv_log_low$x / r1
w1_2 <- inv_log_low$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "inv_log", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "inv_log", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "inv_log", dep_val))

# asym_log cases ------
# high - spherical angle density (using different norms but still spherical angles)
dep_val <- 0.1
asym_log_high <- asym_log(N = 10000, dep = dep_val)
r2 <- sqrt(asym_log_high$x^2 + asym_log_high$y^2)
w2_1 <- asym_log_high$x / r2
w2_2 <- asym_log_high$y / r2
r1 <- asym_log_high$x + asym_log_high$y
w1_1 <- asym_log_high$x / r1
w1_2 <- asym_log_high$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "asym_log", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "asym_log", dep_val, "L1"), col = "green")

# density in terms of L1 norm and L1 angles
hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "asym_log", dep_val))

# asym_log - mid ------
dep_val <- 0.5
asym_log_mid <- asym_log(N = 10000, dep = dep_val)
r2 <- sqrt(asym_log_mid$x^2 + asym_log_mid$y^2)
w2_1 <- asym_log_mid$x / r2
w2_2 <- asym_log_mid$y / r2
r1 <- asym_log_mid$x + asym_log_mid$y
w1_1 <- asym_log_mid$x / r1
w1_2 <- asym_log_mid$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "asym_log", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "asym_log", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "asym_log", dep_val))

# asym_log - low ------
dep_val <- 0.9
asym_log_low <- asym_log(N = 10000, dep = dep_val)
r2 <- sqrt(asym_log_low$x^2 + asym_log_low$y^2)
w2_1 <- asym_log_low$x / r2
w2_2 <- asym_log_low$y / r2
r1 <- asym_log_low$x + asym_log_low$y
w1_1 <- asym_log_low$x / r1
w1_2 <- asym_log_low$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "asym_log", dep_val, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "asym_log", dep_val, "L1"), col = "green")

hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "asym_log", dep_val))

# dirichlet cases ------
# high - spherical angle density (using different norms but still spherical angles)
dep_val <- 1
dep_pair_vals <- c(dep_val, dep_val)
dirichlet_high <- dirichlet(N = 10000, dep_pair_vals)
r2 <- sqrt(dirichlet_high$x^2 + dirichlet_high$y^2)
w2_1 <- dirichlet_high$x / r2
w2_2 <- dirichlet_high$y / r2
r1 <- dirichlet_high$x + dirichlet_high$y
w1_1 <- dirichlet_high$x / r1
w1_2 <- dirichlet_high$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "dirichlet", dep_pair_vals, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "dirichlet", dep_pair_vals, "L1"), col = "green")

# density in terms of L1 norm and L1 angles
hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "dirichlet", dep_pair_vals))

# dirichlet - mid ------
dep_val <- 2
dep_pair_vals <- c(dep_val, dep_val)
dirichlet_high <- dirichlet(N = 10000, dep_pair_vals)
r2 <- sqrt(dirichlet_high$x^2 + dirichlet_high$y^2)
w2_1 <- dirichlet_high$x / r2
w2_2 <- dirichlet_high$y / r2
r1 <- dirichlet_high$x + dirichlet_high$y
w1_1 <- dirichlet_high$x / r1
w1_2 <- dirichlet_high$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "dirichlet", dep_pair_vals, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "dirichlet", dep_pair_vals, "L1"), col = "green")

# density in terms of L1 norm and L1 angles
hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "dirichlet", dep_pair_vals))

# dirichlet - low ------
dep_val <- 0.1
dep_pair_vals <- c(dep_val, dep_val)
dirichlet_high <- dirichlet(N = 10000, dep_pair_vals)
r2 <- sqrt(dirichlet_high$x^2 + dirichlet_high$y^2)
w2_1 <- dirichlet_high$x / r2
w2_2 <- dirichlet_high$y / r2
r1 <- dirichlet_high$x + dirichlet_high$y
w1_1 <- dirichlet_high$x / r1
w1_2 <- dirichlet_high$y / r1
hist(acos(w2_1), freq = FALSE, ylim = c(0,10))
hist(atan(w1_2 / w1_1), freq = FALSE, add = TRUE, col = "orange")
points(acos(w2_1), dens_sphere(w2_1, w2_2, "dirichlet", dep_pair_vals, "L2"), col = "blue")
points(acos(w2_1), dens_sphere(w1_1, w1_2, "dirichlet", dep_pair_vals, "L1"), col = "green")

# density in terms of L1 norm and L1 angles
hist(w1_1, freq = FALSE, ylim = c(0,10))
points(w1_1, dens_l1_norm(w1_1, w1_2, "dirichlet", dep_pair_vals))
