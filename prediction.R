library(geometricMVE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)

## Gauge functions
gauss_gauge <- function(x, y, rho = 0.5) {
  top <- x + y - 2 * rho * sqrt(x * y)
  return(top/(1-rho^2))
}

logistic_gauge <- function(x, y, r = 0.5) {
  r_inv <- 1/r
  return(r_inv * pmax(x, y) + (1-r_inv)*pmin(x,y))
}

inv_log_gauge <- function(x, y, r = 0.5) ((x^(1/r) + y^(1/r))^r)

asym_log_gauge <- function(x, y, r = 0.5) {
  r_inv <- 1/r
  return(pmin((x + y), (r_inv * pmax(x, y) + (1-r_inv)*pmin(x,y))))
}

dirichlet_gauge <- function(x, y, theta1, theta2) {
  return((1 + theta1 + theta2) * pmax(x, y) - (theta1 * x + theta2 * y))
}

rectangular_gauge <- function(x, y, dep) {
  return(pmax((x - y) / dep, (y - x) / dep, (x+ y) / (2 - dep)))
}

# functions ot extact parameters for by gauge function for each dataset
extract_median_params <- function(gauge, dep_type, likelihood, threshold, dep_level, dataset_num) {
  start_file_path <- paste0("stan/csv_fits/stacking/", dep_type, "/", gauge, "/")
  csvfiles <- paste0(start_file_path,
                     list.files(path = start_file_path, 
                                pattern = paste0(dep_level, "_", dataset_num, "_", likelihood, "_", threshold, "_\\d{1}.csv")))
  fit <- as_cmdstan_fit(csvfiles)
  if(gauge != "dirichlet") {
    return(fit |> as_draws_df() |> 
             select(alpha, dep) |> 
             apply(MARGIN = 2, FUN = median) |> t() |>
             as_tibble() |>
             mutate(dataset = dataset_num, gauge = gauge))
  } else {
    return(fit |> as_draws_df() |> 
             select(alpha, theta1, theta2) |> 
             apply(MARGIN = 2, FUN = median) |> t() |>
             as_tibble() |>
             mutate(dataset = dataset_num, gauge = gauge))
  }
}

params_by_gauge_one_dataset <- function(dep_type, likelihood, threshold, dep_level, dataset_num) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  tib_med_params_by_gauge <- lapply(gauge_library, 
                                    function(x) extract_median_params(gauge = x, 
                                                                      dep_type = dep_type,
                                                                      likelihood = likelihood, 
                                                                      threshold = threshold, 
                                                                      dep_level = dep_level, 
                                                                      dataset_num = dataset_num))
  return(tib_med_params_by_gauge)
}

gauss_high_trunc_marg_93 <- params_by_gauge_one_dataset(dep_type = "gauss",
                                                        likelihood = "trunc",
                                                        threshold = "marg",
                                                        dep_level = "high",
                                                        dataset_num = 93)

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


data <- fload("data/gauss/high_93.json")
r0_w <- data$r0_w
X1 <- data$R * data$W
X2 <- data$R * (1-data$W)
R <- data$R
W <- data$W

med_pars <- gauss_high_trunc_marg_93[[1]][,1:2] |> as.numeric()
gw_fit <- gauss_gauge(W, 1-W, med_pars[2])
# gw_true <- gauss_gauge(W, 1-W, 0.9)
plot(W/gw_fit, (1-W)/gw_fit, col = 4)
# points(W/gw_true, (1-W)/gw_true, col = 3)

rw_gw_fit <- cbind(R, W, r0_w, gw_fit) |> as_tibble() |> mutate(r_gw = gw_fit * R) 
ctau <- quantile(rw_gw_fit$r_gw, 0.92)
rw_gw_fit <- rw_gw_fit |> mutate(r0w_tau = ctau / gw_fit)
rw_gw_fit_over1 <- rw_gw_fit |> filter(R > r0w_tau)


r0 <- 4
u <- runif(10000)
q1 <- qgamma(1-u * pgamma(r0, shape = 2, rate = 0.9, lower.tail = FALSE),
             shape = 2, rate = 0.9)
q2 <- qgamma(pgamma(r0, shape = 2, rate = 0.9, lower.tail = TRUE) + u * pgamma(r0, shape = 2, rate = 0.9, lower.tail = FALSE),
             shape = 2, rate = 0.9)
qqplot(q1, q2)

med_pars <- gauss_high_trunc_marg_93[[1]][,1:2] |> as.numeric()

# k = 1
# new_x_liz <- sim.2d_liz(w=rw_gw_fit_over1$W, r0w=rw_gw_fit_over1$r0_w, k=1, 10000, par = med_pars, gfun = gauge_gaussian) |> as_tibble() |> rename(X1 = V1, X2=V2)
new_x <- sim.2d(w=rw_gw_fit$W, r0w=rw_gw_fit$r0w_tau, k=1, 20000, par = med_pars, gfun = gauge_gaussian) |> 
  as_tibble() |> rename(X1 = V1, X2=V2)

plot(X1, X2,pch=20, xlim=c(0,14), ylim = c(0,14))
# points(new_x_liz,pch=20,col=4)
points(new_x,pch=20,col=3)

lower <- 8
upper <- 10
# prob_new_x_liz <- (new_x_liz |> filter(X1 >= lower & X1 <= upper & X2 >= lower & X2 <= upper) |> nrow())/nrow(new_x_liz) * nrow(rw_gw_fit_over1)/nrow(rw_gw_fit)
(new_x |> filter(X1 >= lower & X1 <= upper & X2 >= lower & X2 <= upper) |> 
                 nrow())/nrow(new_x) * 0.08

# convert [8,10] x [8,10] box in expo to gaussian coordinates
lower_limit <- qnorm(pexp(lower))
upper_limit <- qnorm(pexp(upper))
pmvnorm(lower = rep(lower_limit, 2), upper = rep(upper_limit, 2), corr = matrix(c(1, 0.9, 0.9, 1), nrow = 2))[1]


