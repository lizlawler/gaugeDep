library(gaugeDependence)

data <- RcppSimdJson::fload("../gaugeDep/data/gauss/low_20.json")
w <- data$W
r <- data$R
r0w <- data$r0_w
idx <- data$idx



temp <- seq(0, 1, length.out = 100)
grid <- expand.grid(temp, temp)
# sum_grid <- grid[,1] + grid[,2]
# sqrt_grid <- sqrt(grid[,1] * grid[,2])

angle_results <- angular_mcmc(angles = w, dim = 2, 
                              starting_theta = runif(1),
                              # starting_theta = c(abs(rt(1, 4,ncp = 0))*4, abs(rt(1, 4,ncp = 0))*2), 
                              gauge_type = "inv_log", 
                              n_updates = 15000, 
                              update_freq = 250, 
                              n_burnin = 5000,
                              n_thin = 5,
                              adapt_cov = TRUE)
trunc_results <- radial_adaptive_mh(radii = r[idx], r0w = r0w[idx], angles = w[idx],
                              starting_theta = c(rgamma(1, 4, 2), runif(1)),
                              # starting_theta = c(rgamma(1, 4, 2), abs(rt(1, 4,ncp = 0)), abs(rt(1, 4,ncp = 0))),
                              likelihood_type = "trunc",
                              gauge_type = "gauss",
                              n_updates = 15000, 
                              update_freq = 250, 
                              n_burnin = 5000, 
                              n_thin = 5, 
                              adapt_cov = TRUE)

cens_results <- radial_adaptive_mh(radii = r, r0w = r0w, angles = w,
                                   starting_theta = c(rgamma(1, 4, 2), runif(1)),
                                   # starting_theta = c(rgamma(1, 4, 2), abs(rt(1, 4,ncp = 0)), abs(rt(1, 4,ncp = 0))),
                                   likelihood_type = "cens",
                                   gauge_type = "gauss",
                                   n_updates = 15000, 
                                   update_freq = 250, 
                                   n_burnin = 5000,
                                   n_thin = 5,
                                   adapt_cov = TRUE)
plot(cens_results$sample[,"dep"], type = "l")
plot(density(cens_results$samples[,"dep"]))
plot(angle_results$sample[,"dep"], type = "l")
lines(cens_results$sample[,"dep"], col = "blue")

w <- seq(0, 1, length.out = 200)
gw <- gaugeDependence::dirichlet_gauge(w, c(3, 1))
plot(w/gw, (1-w)/gw)
