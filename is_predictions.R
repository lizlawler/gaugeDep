library(tidyverse)
library(cmdstanr)
source("gauge_functions_wrt_w.R")
source("gauge_functions_wrt_x.R")

trunc_gamma <- function(x, xmin, alpha, beta) {
  unnorm_pdf <- dgamma(x, shape = alpha, rate = beta)
  norm_cst <- pgamma(xmin, shape = alpha, rate = beta, lower.tail = F)
  return(unnorm_pdf / norm_cst)
}

est_volume <- function(n = 100, pars = 0.5, gauge) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  gauge_fcn <- get(paste0(gauge, "_gauge_wrt_x"))
  gx <- gauge_fcn(grid[,1], grid[,2], dep_par = pars)
  return(mean(gx <= 1))
}

dens_l1_norm <- function(w1, gauge, par_val) {
  if(gauge == "logistic") {
    mc_vol <- par_val
  } else {
    mc_vol <- est_volume(n = 100, par_val, gauge)
  }
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  gw <- gauge_fcn(w1, par_val)
  return(1 / (gw^2 * 2 * mc_vol))
}

true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(mvtnorm::pmvnorm(lower = c(dim1_star[1],dim2_star[1]), upper = c(dim1_star[2],dim2_star[2]), corr = corr_matrix)[1])
}

#### HIGH DEPENDENCE ####
data <- RcppSimdJson::fload("data/gauss/high_1.json")
r <- data$R
w <- data$W
r0w <- data$r0_w_ctau

gauss_fit <- readRDS("~/Desktop/research/gaugeDependence/mcmc_samples/gauss_high_1.RDS")
post_alpha <- median(gauss_fit$trace[10000:50000, 1])
post_theta1 <- median(gauss_fit$trace[10000:50000, 2])
post_theta2 <- median(gauss_fit$trace[10000:50000, 3])
# IS sampling for R|W * W - method 1 (using MVN)
is_samp_mvn <- mvtnorm::rmvnorm(5000, c(11,11), sigma = 2 * diag(2)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_mvn <- is_samp_mvn |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_mvn$w1, post_theta2))
is_samp_rw_mvn <- is_samp_rw_mvn |> mutate(r0w = r0w_post)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_mvn, pch = 20, col = alpha("red", 0.5))

w_subset_mvn <- is_samp_rw_mvn |> filter(w1 >= 5/11, w1 <= 6/11)
points(w_subset_mvn[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 5/6)
abline(a = 0, b = 6/5)
abline(v = 10)
abline(v = 12)
abline(h = 12)
abline(h = 10)
r_subset_w_mvn <- w_subset_mvn |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (10/(1-w1)), r <= (12/(1-w1)))
points(r_subset_w_mvn[,1:2], pch = 20, col = "green")

# probability calculation
rw_dens_mvn <- trunc_gamma(r_subset_w_mvn$r, r_subset_w_mvn$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_mvn$w1, post_theta2)) *
  dens_l1_norm(r_subset_w_mvn$w1, "gauss", post_theta1)
is_dens_mvn <- mvtnorm::dmvnorm(r_subset_w_mvn[,1:2], mean = c(11,11), sigma = 2 * diag(2)) * r_subset_w_mvn$r
weights_mvn <- rw_dens_mvn / is_dens_mvn
(sum(weights_mvn) / nrow(is_samp_rw_mvn)) * 0.05

true_gauss_prob(c(10,12), c(10,12), 0.5)

# IS sampling for R|W, in B1 - method 2 (using Gamma)
is_samp_gamma <- cbind(rgamma(5000, 11, 1), rgamma(5000, 11, 1)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_gamma <- is_samp_gamma |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_gamma$w1, post_theta2))
is_samp_rw_gamma <- is_samp_rw_gamma |> mutate(r0w = r0w_post)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_gamma, pch = 20, col = alpha("red", 0.5))

w_subset_gamma <- is_samp_rw_gamma |> filter(w1 >= 5/11, w1 <= 6/11)
points(w_subset_gamma[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 5/6)
abline(a = 0, b = 6/5)
abline(v = 10)
abline(v = 12)
abline(h = 12)
abline(h = 10)
r_subset_w_gamma <- w_subset_gamma |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (10/(1-w1)), r <= (12/(1-w1)))
points(r_subset_w_gamma[,1:2], pch = 20, col = "green")

r_giv_w_dens_gamma <- trunc_gamma(r_subset_w_gamma$r, r_subset_w_gamma$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_gamma$w1, post_theta2))
is_dens_cond_gamma <- dgamma(r_subset_w_gamma$x1, 11, 1) * dgamma(r_subset_w_gamma$x2, 11, 1) * r_subset_w_gamma$r / dbeta(r_subset_w_gamma$w1, 11, 11)
weights_r <- r_giv_w_dens_gamma / is_dens_cond_gamma
prob_r_box <- sum(weights_r) / nrow(is_samp_gamma)

# IS sampling for W in B1
w_is <- runif(1000)
is_idx <- which(w_is <= 6/11 & w_is >= 5/11)
weights_w <- dens_l1_norm(w_is[is_idx], "gauss", par_val = post_theta1)
prob_w_box <- sum(weights_w)/1000

prob_r_box * prob_w_box * 0.05

# IS sampling for R|W * W - method 1 (using MVN), box 2
is_samp_mvn <- mvtnorm::rmvnorm(5000, c(11,7), sigma = 2 * diag(2)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_mvn <- is_samp_mvn |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_mvn$w1, post_theta2))
is_samp_rw_mvn <- is_samp_rw_mvn |> mutate(r0w = r0w_post)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_mvn, pch = 20, col = alpha("red", 0.5))

w_subset_mvn <- is_samp_rw_mvn |> filter(w1 >= 5/9, w1 <= 6/9)
points(w_subset_mvn[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 4/5)
abline(a = 0, b = 1/2)
abline(v = 10)
abline(v = 12)
abline(h = 8)
abline(h = 6)
r_subset_w_mvn <- w_subset_mvn |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (6/(1-w1)), r <= (8/(1-w1)))
points(r_subset_w_mvn[,1:2], pch = 20, col = "green")

# probability calculation
rw_dens_mvn <- trunc_gamma(r_subset_w_mvn$r, r_subset_w_mvn$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_mvn$w1, post_theta2)) *
  dens_l1_norm(r_subset_w_mvn$w1, "gauss", post_theta1)
is_dens_mvn <- mvtnorm::dmvnorm(r_subset_w_mvn[,1:2], mean = c(11,7), sigma = 2 * diag(2)) * r_subset_w_mvn$r
weights_mvn <- rw_dens_mvn / is_dens_mvn
(sum(weights_mvn) / nrow(is_samp_rw_mvn)) * 0.05

true_gauss_prob(c(10,12), c(6,8), 0.5)

# IS sampling for R|W, in B1 - method 2 (using Gamma)
is_samp_gamma <- cbind(rgamma(5000, 12, 1), rgamma(5000, 8, 1)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_gamma <- is_samp_gamma |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_gamma$w1, post_theta2))
is_samp_rw_gamma <- is_samp_rw_gamma |> mutate(r0w = r0w_post)

# check all points are above threshold
sum(is_samp_rw_gamma$r > is_samp_rw_gamma$r0w)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_gamma, pch = 20, col = alpha("red", 0.5))

w_subset_gamma <- is_samp_rw_gamma |> filter(w1 >= 5/9, w1 <= 6/9)
points(w_subset_gamma[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 4/5)
abline(a = 0, b = 1/2)
abline(v = 10)
abline(v = 12)
abline(h = 6)
abline(h = 8)
r_subset_w_gamma <- w_subset_gamma |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (6/(1-w1)), r <= (8/(1-w1)))
points(r_subset_w_gamma[,1:2], pch = 20, col = "green")

r_giv_w_dens_gamma <- trunc_gamma(r_subset_w_gamma$r, r_subset_w_gamma$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_gamma$w1, post_theta2))
is_dens_cond_gamma <- dgamma(r_subset_w_gamma$x1, 12, 1) * dgamma(r_subset_w_gamma$x2, 8, 1) * r_subset_w_gamma$r / dbeta(r_subset_w_gamma$w1, 12, 8)
weights_r <- r_giv_w_dens_gamma / is_dens_cond_gamma
prob_r_box <- sum(weights_r) / nrow(is_samp_gamma)

# IS sampling for W in B2
w_is <- runif(1000)
is_idx <- which(w_is <= 6/9 & w_is >= 5/9)
weights_w <- dens_l1_norm(w_is[is_idx], "gauss", par_val = post_theta1)
prob_w_box <- sum(weights_w)/1000

prob_r_box * prob_w_box * 0.05
true_gauss_prob(c(10,12), c(6,8), 0.5)

# IS sampling for R|W * W - method 1 (using MVN), box 3
is_samp_mvn <- mvtnorm::rmvnorm(5000, c(11,3), sigma = 2 * diag(2)) |> as_tibble() |> rename(x1 = V1, x2 = V2) |> filter(x1 >= 0, x2 >= 0)
is_samp_rw_mvn <- is_samp_mvn |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_mvn$w1, post_theta2))
is_samp_rw_mvn <- is_samp_rw_mvn |> mutate(r0w = r0w_post)

sum(is_samp_rw_mvn$r > is_samp_rw_mvn$r0w)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_mvn, pch = 20, col = alpha("red", 0.5))

w_subset_mvn <- is_samp_rw_mvn |> filter(w1 >= 5/7, w1 <= 6/7)
points(w_subset_mvn[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 2/5)
abline(a = 0, b = 1/6)
abline(v = 10)
abline(v = 12)
abline(h = 4)
abline(h = 2)
r_subset_w_mvn <- w_subset_mvn |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (2/(1-w1)), r <= (4/(1-w1)))
points(r_subset_w_mvn[,1:2], pch = 20, col = "green")

# probability calculation
rw_dens_mvn <- trunc_gamma(r_subset_w_mvn$r, r_subset_w_mvn$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_mvn$w1, post_theta2)) *
  dens_l1_norm(r_subset_w_mvn$w1, "gauss", post_theta1)
is_dens_mvn <- mvtnorm::dmvnorm(r_subset_w_mvn[,1:2], mean = c(11,3), sigma = 2 * diag(2)) * r_subset_w_mvn$r
weights_mvn <- rw_dens_mvn / is_dens_mvn
(sum(weights_mvn) / nrow(is_samp_rw_mvn)) * 0.05

true_gauss_prob(c(10,12), c(2,4), 0.5)

# IS sampling for R|W, in B1 - method 3 (using Gamma)
is_samp_gamma <- cbind(rgamma(5000, 12, 1), rgamma(5000, 4, 1)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_gamma <- is_samp_gamma |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_gamma$w1, post_theta2))
is_samp_rw_gamma <- is_samp_rw_gamma |> mutate(r0w = r0w_post)

# check all points are above threshold
sum(is_samp_rw_gamma$r > is_samp_rw_gamma$r0w)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_gamma, pch = 20, col = alpha("red", 0.5))

w_subset_gamma <- is_samp_rw_gamma |> filter(w1 >= 5/7, w1 <= 6/7)
points(w_subset_gamma[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 2/5)
abline(a = 0, b = 1/6)
abline(v = 10)
abline(v = 12)
abline(h = 4)
abline(h = 2)
r_subset_w_gamma <- w_subset_gamma |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (2/(1-w1)), r <= (4/(1-w1)))
points(r_subset_w_gamma[,1:2], pch = 20, col = "green")

r_giv_w_dens_gamma <- trunc_gamma(r_subset_w_gamma$r, r_subset_w_gamma$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_gamma$w1, post_theta2))
is_dens_cond_gamma <- dgamma(r_subset_w_gamma$x1, 12, 1) * dgamma(r_subset_w_gamma$x2, 4, 1) * r_subset_w_gamma$r / dbeta(r_subset_w_gamma$w1, 12, 4)
weights_r <- r_giv_w_dens_gamma / is_dens_cond_gamma
prob_r_box <- sum(weights_r) / nrow(is_samp_gamma)

# IS sampling for W in B1
w_is <- runif(1000)
is_idx <- which(w_is <= 6/7 & w_is >= 5/7)
weights_w <- dens_l1_norm(w_is[is_idx], "gauss", par_val = post_theta1)
prob_w_box <- sum(weights_w)/1000

prob_r_box * prob_w_box * 0.05
true_gauss_prob(c(10,12), c(2,4), 0.5)


#### MID DEPENDENCE ####
data <- RcppSimdJson::fload("data/gauss/mid_1.json")
r <- data$R
w <- data$W
r0w <- data$r0_w_ctau

gauss_fit <- readRDS("~/Desktop/research/gaugeDependence/mcmc_samples/gauss_high_1.RDS")
post_alpha <- median(gauss_fit$trace[10000:50000, 1])
post_theta1 <- median(gauss_fit$trace[10000:50000, 2])
post_theta2 <- median(gauss_fit$trace[10000:50000, 3])
# IS sampling for R|W * W - method 1 (using MVN)
is_samp_mvn <- mvtnorm::rmvnorm(5000, c(11,11), sigma = 2 * diag(2)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_mvn <- is_samp_mvn |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_mvn$w1, post_theta2))
is_samp_rw_mvn <- is_samp_rw_mvn |> mutate(r0w = r0w_post)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_mvn, pch = 20, col = alpha("red", 0.5))

w_subset_mvn <- is_samp_rw_mvn |> filter(w1 >= 5/11, w1 <= 6/11)
points(w_subset_mvn[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 5/6)
abline(a = 0, b = 6/5)
abline(v = 10)
abline(v = 12)
abline(h = 12)
abline(h = 10)
r_subset_w_mvn <- w_subset_mvn |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (10/(1-w1)), r <= (12/(1-w1)))
points(r_subset_w_mvn[,1:2], pch = 20, col = "green")

# probability calculation
rw_dens_mvn <- trunc_gamma(r_subset_w_mvn$r, r_subset_w_mvn$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_mvn$w1, post_theta2)) *
  dens_l1_norm(r_subset_w_mvn$w1, "gauss", post_theta1)
is_dens_mvn <- mvtnorm::dmvnorm(r_subset_w_mvn[,1:2], mean = c(11,11), sigma = 2 * diag(2)) * r_subset_w_mvn$r
weights_mvn <- rw_dens_mvn / is_dens_mvn
(sum(weights_mvn) / nrow(is_samp_rw_mvn)) * 0.05

true_gauss_prob(c(10,12), c(10,12), 0.5)

# IS sampling for R|W, in B1 - method 2 (using Gamma)
is_samp_gamma <- cbind(rgamma(5000, 11, 1), rgamma(5000, 11, 1)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_gamma <- is_samp_gamma |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_gamma$w1, post_theta2))
is_samp_rw_gamma <- is_samp_rw_gamma |> mutate(r0w = r0w_post)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_gamma, pch = 20, col = alpha("red", 0.5))

w_subset_gamma <- is_samp_rw_gamma |> filter(w1 >= 5/11, w1 <= 6/11)
points(w_subset_gamma[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 5/6)
abline(a = 0, b = 6/5)
abline(v = 10)
abline(v = 12)
abline(h = 12)
abline(h = 10)
r_subset_w_gamma <- w_subset_gamma |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (10/(1-w1)), r <= (12/(1-w1)))
points(r_subset_w_gamma[,1:2], pch = 20, col = "green")

r_giv_w_dens_gamma <- trunc_gamma(r_subset_w_gamma$r, r_subset_w_gamma$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_gamma$w1, post_theta2))
is_dens_cond_gamma <- dgamma(r_subset_w_gamma$x1, 11, 1) * dgamma(r_subset_w_gamma$x2, 11, 1) * r_subset_w_gamma$r / dbeta(r_subset_w_gamma$w1, 11, 11)
weights_r <- r_giv_w_dens_gamma / is_dens_cond_gamma
prob_r_box <- sum(weights_r) / nrow(is_samp_gamma)

# IS sampling for W in B1
w_is <- runif(1000)
is_idx <- which(w_is <= 6/11 & w_is >= 5/11)
weights_w <- dens_l1_norm(w_is[is_idx], "gauss", par_val = post_theta1)
prob_w_box <- sum(weights_w)/1000

prob_r_box * prob_w_box * 0.05
true_gauss_prob(c(10,12), c(10,12), 0.5)

# IS sampling for R|W * W - method 1 (using MVN), box 2
is_samp_mvn <- mvtnorm::rmvnorm(5000, c(11,7), sigma = 2 * diag(2)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_mvn <- is_samp_mvn |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_mvn$w1, post_theta2))
is_samp_rw_mvn <- is_samp_rw_mvn |> mutate(r0w = r0w_post)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_mvn, pch = 20, col = alpha("red", 0.5))

w_subset_mvn <- is_samp_rw_mvn |> filter(w1 >= 5/9, w1 <= 6/9)
points(w_subset_mvn[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 4/5)
abline(a = 0, b = 1/2)
abline(v = 10)
abline(v = 12)
abline(h = 8)
abline(h = 6)
r_subset_w_mvn <- w_subset_mvn |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (6/(1-w1)), r <= (8/(1-w1)))
points(r_subset_w_mvn[,1:2], pch = 20, col = "green")

# probability calculation
rw_dens_mvn <- trunc_gamma(r_subset_w_mvn$r, r_subset_w_mvn$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_mvn$w1, post_theta2)) *
  dens_l1_norm(r_subset_w_mvn$w1, "gauss", post_theta1)
is_dens_mvn <- mvtnorm::dmvnorm(r_subset_w_mvn[,1:2], mean = c(11,7), sigma = 2 * diag(2)) * r_subset_w_mvn$r
weights_mvn <- rw_dens_mvn / is_dens_mvn
(sum(weights_mvn) / nrow(is_samp_rw_mvn)) * 0.05

true_gauss_prob(c(10,12), c(6,8), 0.5)

# IS sampling for R|W, in B1 - method 2 (using Gamma)
is_samp_gamma <- cbind(rgamma(5000, 12, 1), rgamma(5000, 8, 1)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_gamma <- is_samp_gamma |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_gamma$w1, post_theta2))
is_samp_rw_gamma <- is_samp_rw_gamma |> mutate(r0w = r0w_post)

# check all points are above threshold
sum(is_samp_rw_gamma$r > is_samp_rw_gamma$r0w)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_gamma, pch = 20, col = alpha("red", 0.5))

w_subset_gamma <- is_samp_rw_gamma |> filter(w1 >= 5/9, w1 <= 6/9)
points(w_subset_gamma[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 4/5)
abline(a = 0, b = 1/2)
abline(v = 10)
abline(v = 12)
abline(h = 6)
abline(h = 8)
r_subset_w_gamma <- w_subset_gamma |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (6/(1-w1)), r <= (8/(1-w1)))
points(r_subset_w_gamma[,1:2], pch = 20, col = "green")

r_giv_w_dens_gamma <- trunc_gamma(r_subset_w_gamma$r, r_subset_w_gamma$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_gamma$w1, post_theta2))
is_dens_cond_gamma <- dgamma(r_subset_w_gamma$x1, 12, 1) * dgamma(r_subset_w_gamma$x2, 8, 1) * r_subset_w_gamma$r / dbeta(r_subset_w_gamma$w1, 12, 8)
weights_r <- r_giv_w_dens_gamma / is_dens_cond_gamma
prob_r_box <- sum(weights_r) / nrow(is_samp_gamma)

# IS sampling for W in B2
w_is <- runif(1000)
is_idx <- which(w_is <= 6/9 & w_is >= 5/9)
weights_w <- dens_l1_norm(w_is[is_idx], "gauss", par_val = post_theta1)
prob_w_box <- sum(weights_w)/1000

prob_r_box * prob_w_box * 0.05
true_gauss_prob(c(10,12), c(6,8), 0.5)

# IS sampling for R|W * W - method 1 (using MVN), box 3
is_samp_mvn <- mvtnorm::rmvnorm(5000, c(11,3), sigma = 2 * diag(2)) |> as_tibble() |> rename(x1 = V1, x2 = V2) |> filter(x1 >= 0, x2 >= 0)
is_samp_rw_mvn <- is_samp_mvn |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_mvn$w1, post_theta2))
is_samp_rw_mvn <- is_samp_rw_mvn |> mutate(r0w = r0w_post)

sum(is_samp_rw_mvn$r > is_samp_rw_mvn$r0w) / nrow(is_samp_rw_mvn)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_mvn, pch = 20, col = alpha("red", 0.5))

w_subset_mvn <- is_samp_rw_mvn |> filter(w1 >= 5/7, w1 <= 6/7)
points(w_subset_mvn[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 2/5)
abline(a = 0, b = 1/6)
abline(v = 10)
abline(v = 12)
abline(h = 4)
abline(h = 2)
r_subset_w_mvn <- w_subset_mvn |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (2/(1-w1)), r <= (4/(1-w1)))
points(r_subset_w_mvn[,1:2], pch = 20, col = "green")

# probability calculation
rw_dens_mvn <- trunc_gamma(r_subset_w_mvn$r, r_subset_w_mvn$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_mvn$w1, post_theta2)) *
  dens_l1_norm(r_subset_w_mvn$w1, "gauss", post_theta1)
is_dens_mvn <- mvtnorm::dmvnorm(r_subset_w_mvn[,1:2], mean = c(11,3), sigma = 2 * diag(2)) * r_subset_w_mvn$r
weights_mvn <- rw_dens_mvn / is_dens_mvn
(sum(weights_mvn) / nrow(is_samp_rw_mvn)) * 0.05

true_gauss_prob(c(10,12), c(2,4), 0.5)

# IS sampling for R|W, in B1 - method 3 (using Gamma)
is_samp_gamma <- cbind(rgamma(5000, 12, 1), rgamma(5000, 4, 1)) |> as_tibble() |> rename(x1 = V1, x2 = V2)
is_samp_rw_gamma <- is_samp_gamma |> mutate(r = x1 + x2, w1 = x1/r)
r0w_post <- qgamma(0.95, shape = post_alpha, rate = gauss_gauge(is_samp_rw_gamma$w1, post_theta2))
is_samp_rw_gamma <- is_samp_rw_gamma |> mutate(r0w = r0w_post)

# check all points are above threshold
sum(is_samp_rw_gamma$r > is_samp_rw_gamma$r0w)

plot(w * r, (1-w) * r, pch = 20, xlim = c(0,14), ylim = c(0,14))
points(is_samp_gamma, pch = 20, col = alpha("red", 0.5))

w_subset_gamma <- is_samp_rw_gamma |> filter(w1 >= 5/7, w1 <= 6/7)
points(w_subset_gamma[,1:2], pch = 20, col = "blue")
abline(a = 0, b = 2/5)
abline(a = 0, b = 1/6)
abline(v = 10)
abline(v = 12)
abline(h = 4)
abline(h = 2)
r_subset_w_gamma <- w_subset_gamma |> 
  filter(r >= (10/w1), r <= (12/w1), r >= (2/(1-w1)), r <= (4/(1-w1)))
points(r_subset_w_gamma[,1:2], pch = 20, col = "green")

r_giv_w_dens_gamma <- trunc_gamma(r_subset_w_gamma$r, r_subset_w_gamma$r0w, alpha = post_alpha, beta = gauss_gauge(r_subset_w_gamma$w1, post_theta2))
is_dens_cond_gamma <- dgamma(r_subset_w_gamma$x1, 12, 1) * dgamma(r_subset_w_gamma$x2, 4, 1) * r_subset_w_gamma$r / dbeta(r_subset_w_gamma$w1, 12, 4)
weights_r <- r_giv_w_dens_gamma / is_dens_cond_gamma
prob_r_box <- sum(weights_r) / nrow(is_samp_gamma)

# IS sampling for W in B1
w_is <- runif(1000)
is_idx <- which(w_is <= 6/7 & w_is >= 5/7)
weights_w <- dens_l1_norm(w_is[is_idx], "gauss", par_val = post_theta1)
prob_w_box <- sum(weights_w)/1000

prob_r_box * prob_w_box * 0.05
true_gauss_prob(c(10,12), c(2,4), 0.5)

rw_dens_gamma <- trunc_gamma(r_subset_w_gamma$r, r_subset_w_gamma$r0w, 
                             alpha = post_alpha, beta = gauss_gauge(r_subset_w_gamma$w1, post_theta2)) *
  dens_l1_norm(r_subset_w_gamma$w1, gauge = "gauss", post_theta1)
is_dens_gamma <- dgamma(r_subset_w_gamma$x1, 12, 1) * dgamma(r_subset_w_gamma$x2, 4, 1) * r_subset_w_gamma$r
weights_rw_gamma <- rw_dens_gamma / is_dens_gamma
(sum(weights_rw_gamma) / nrow(is_samp_rw_gamma)) * 0.05

