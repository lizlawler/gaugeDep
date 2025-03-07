library(BezELS)
library(RcppSimdJson)
library(tidyverse)
library(gaugeDependence)

data <- fload("data/gauss/high_1.json")
idx <- data$idx
samples  <- fit_mcmc_bezier( N = data$n0,
                             r = data$R[idx],
                             w = data$W[idx],
                             r_0 = data$r0_w[idx],
                             iters = 11000, burn = 100, 
                             traceplot = FALSE,
                             print.every = 100)

x <- cbind(data$R * data$W, data$R * (1 - data$W))
r0w <- data$r0_w
above <- data$R > r0w

plot_bezier_polar(samples, x = x, r_0_marg = r0w, above_thresh_marg = above, theta = 0.9, copula = "g")

angles <- cbind(data$W[idx], 1 - data$W[idx])

post_gamma_rate <- function(mcmc_samples, w_mat, thin.by = 10) {
  n.data = nrow(w_mat)
  n.samples = nrow(mcmc_samples)
  sss = seq(1, n.samples, by = thin.by)
  
  n.post = length(sss)
  post_rate <- matrix(NA, nrow = n.post, ncol = n.data)
  
  p0y <- (mcmc_samples[sss, 1])
  p1x <- (mcmc_samples[sss, 2])
  p1y <- (mcmc_samples[sss, 3])
  p2x <- (mcmc_samples[sss, 4])
  p3 <- (mcmc_samples[sss, 5])
  p4y <- (mcmc_samples[sss, 6])
  p5x <- (mcmc_samples[sss, 7])
  p5y <- (mcmc_samples[sss, 8])
  p6x <- (mcmc_samples[sss, 9])
  
  for (i in 1:n.post) {
    p <- matrix(c(0, p0y[i], p1x[i], p1y[i], p2x[i], 1, p3[i], 
                  p3[i], 1, p4y[i], p5x[i], p5y[i], p6x[i], 0), 
                nrow = 7, ncol = 2, byrow = TRUE)
    post_rate[i, ] <- gx(N = n.data, p, w_mat[, 2]/w_mat[, 1], w_mat[, 1])
  }
  
  return(apply(post_rate, 2, median))
}

gamma_rate <- post_gamma_rate(samples, angles, 10)

test <- qs::qread("fits_and_weights/post_params_joint/gauss_gauss_high_cens_radial.qs")
gamma_rate_nobez <- as.numeric(gauss_gauge(angles[,1], angles[,2], test$dep[1]))

plot(angles[,1], gamma_rate, pch = 20, col = "blue")
points(angles[,1], gamma_rate_nobez, pch = 20, col = "red")

data <- fload("data/gauss/low_1.json")
idx <- data$idx
samples  <- fit_mcmc_bezier( N = data$n0,
                             r = data$R[idx],
                             w = data$W[idx],
                             r_0 = data$r0_w[idx],
                             iters = 11000, burn = 100, 
                             traceplot = FALSE,
                             print.every = 100)

angles <- cbind(data$W[idx], 1 - data$W[idx])
gamma_rate <- post_gamma_rate(samples, angles, 10)

test <- qs::qread("fits_and_weights/post_params_joint/gauss_gauss_low_cens_radial.qs")
gamma_rate_nobez <- as.numeric(gauss_gauge(angles[,1], angles[,2], test$dep[1]))

plot(angles[,1], gamma_rate, pch = 20, col = "blue")
points(angles[,1], gamma_rate_nobez, pch = 20, col = "red")

x <- cbind(data$R * data$W, data$R * (1 - data$W))
r0w <- data$r0_w
above <- data$R > r0w

plot_bezier_polar(samples, x = x, r_0_marg = r0w, above_thresh_marg = above, theta = 0.1, copula = "g")
