# library(geometricMVE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)

## Gauge functions
gauss_gauge <- function(w, dep_par = 0.5) {
  top <- w + (1 - w) - 2 * dep_par * sqrt(w * (1 - w))
  return(top/(1-dep_par^2))
}

logistic_gauge <- function(w, dep_par = 0.5) {
  r_inv <- 1/dep_par
  return(r_inv * pmax(w, (1 - w)) + (1-r_inv)*pmin(w,(1 - w)))
}

inv_log_gauge <- function(w, dep_par = 0.5) ((w^(1/dep_par) + (1 - w)^(1/dep_par))^dep_par)

asym_log_gauge <- function(w, dep_par = 0.5) {
  r_inv <- 1/dep_par
  return(pmin((w + (1 - w)), (r_inv * pmax(w, (1 - w)) + (1-r_inv)*pmin(w,(1 - w)))))
}

dirichlet_gauge <- function(w, dep_par) {
  theta1 <- dep_par[1]
  theta2 <- dep_par[2]
  return((1 + theta1 + theta2) * pmax(w, (1 - w)) - (theta1 * w + theta2 * (1 - w)))
}

rectangular_gauge <- function(w, dep_par) {
  return(pmax((w - (1 - w)) / dep_par, ((1 - w) - w) / dep_par, (w + (1 - w)) / (2 - dep_par)))
}

# create function to simulate new data when k = 1 -----------
sim_new_data <- function(w, r0w, nsim, pars, gauge_fcn) {
  resample_idx <- sample(1:length(w), size = nsim, replace = T)
  w_resample <- w[resample_idx]
  r0w_resample <- r0w[resample_idx]
  rate_init <- gauge_fcn(w, dep_par = pars[2:length(pars)])
  rate_star <- rate_init[resample_idx] 
  rstar <- qgamma(1 - runif(nsim) * pgamma(r0w_resample, shape = pars[1], 
                                           rate = rate_star, lower.tail = F),
                  shape = pars[1], rate = rate_star)
  xstar <- cbind(rstar * w_resample, rstar * (1 - w_resample))
  return(xstar |> as_tibble() |> rename(x1 = V1, x2 = V2))
}

# create function to make predictions from newly simulated data -------
pred_probs <- function(sim_df, lower, upper, idx, length_data) {
  prob_x_given_r <- (sim_df |>
                       filter(x1 >= lower & x1 <= upper & x2 >= lower & x2 <= upper) |>
                       nrow()) / nrow(sim_df)
  prob_r_over_thres <- length(idx)/length_data
  return(prob_x_given_r * prob_r_over_thres)
}

# create function that wraps everything together ---------
make_preds <- function(data_file, lower, upper, posterior_pars, gauge) {
  data <- fload(data_file)
  r0_w <- data$r0_w
  X1 <- data$R * data$W
  X2 <- data$R * (1-data$W)
  R <- data$R
  W <- data$W
  idx <- data$idx
  sim_df <- sim_new_data(w = W[idx], r0w = r0_w[idx], nsim = 10000, pars = posterior_pars, gauge_fcn = get(paste0(gauge, "_gauge")))
  return(pred_probs(sim_df, lower, upper, idx, length(data$R)))
}

# create function to make predictions by the gauge function it was fit to ----------
preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, threshold) {
  if(dep_level == "high") {
    lower_lim <- 8
    upper_lim <- 10
  } else if (dep_level == "mid") {
    lower_lim <- 7
    upper_lim <- 9
  } else {
    lower_lim <- 6
    upper_lim <- 8
  }
  posterior_params <- readRDS(paste0("extracted_params/", gauge, "_", dep_type, "_", 
                                     dep_level, "_", likelihood, "_", threshold, "_params.RDS"))
  init_data_path <- paste0("data/", dep_type, "/", dep_level, "_")
  results <- apply(posterior_params, 1, function(row) {
    data_file <- paste0(init_data_path, row["dataset"], ".json")
    params <- as.numeric(row[-length(row)])
    return(make_preds(data_file, lower = lower_lim, upper = upper_lim, posterior_pars = params, gauge = gauge))
  })
  return(results)
}

# create function that makes predictions for all 100 datasets for a specific dependence type and level, likelihood type,
# threshold type, and with all gauge function fits -------
preds_by_dep_level_lhood_thres <- function(dep_type, dep_level, likelihood, threshold) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(sapply(gauge_library, function(x) preds_by_gauge(x, dep_type, dep_level, likelihood, threshold)) |>
           as_tibble() |>
           mutate(dataset = 1:100))
}

# create function to reshape previously extracted stacking weights -------
make_wts_df <- function(weights_file) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  temp <- readRDS(weights_file) |>
    bind_rows() |> 
    mutate(method = rep(gauge_library, 100)) |>
    mutate(stacking = as.numeric(stacking),
           pseudobma_boot = as.numeric(pseudobma_boot),
           pseudobma_noboot = as.numeric(pseudobma_noboot)) |>
    mutate(dataset = rep(1:100, times = rep(6, 100)))
  return(temp)
}

# create function to create weighted sum of predictions by three BMA methods ------
weighted_preds_by_lhood_thres <- function(dep_type, dep_level, likelihood, threshold) {
  scenario <- paste0(dep_type, "_", dep_level, "_", likelihood, "_", threshold)
  temp_preds <- preds_by_dep_level_lhood_thres(dep_type, dep_level, likelihood, threshold) |>
    pivot_longer(cols = -'dataset', names_to = "method", values_to = "preds")
  temp_wts <- make_wts_df(paste0("stacking_weights/", scenario, "_wts.RDS"))
  temp_weighted_preds <- temp_wts |> left_join(temp_preds) |>
    mutate(stacking_preds = preds * stacking,
           pseudo_boot = pseudobma_boot * preds,
           pseudo_noboot = pseudobma_noboot * preds) |>
    group_by(dataset) |>
    summarize(stacking_predictions = sum(stacking_preds),
              pseudobma_boot_preds = sum(pseudo_boot),
              pseudobma_noboot_preds = sum(pseudo_noboot)) |>
    ungroup()
  boxplot_wts <- temp_weighted_preds |> 
    pivot_longer(cols = -'dataset', names_to = "bma_method", values_to = "bma_preds") |>
    mutate(bma_method = case_when(grepl("stacking", bma_method) ~ 'Stacking',
                                  grepl("noboot", bma_method) ~ 'Pseudo-BMA',
                                  grepl("boot", bma_method) ~ 'Pseudo-BMA+'),
           bma_method = as.factor(bma_method),
           scenario = rep(scenario))
  return(boxplot_wts)
}

weighted_preds_by_level <- function(dep_type, dep_level) {
  thres <- c("ctau", "marg")
  lhood <- c("trunc", "cens")
  lhood_thres_combos <- expand_grid(lhood, thres)
  all_wts <- apply(lhood_thres_combos, 1, 
                   function(row) weighted_preds_by_lhood_thres(dep_type, dep_level, 
                                                               row["lhood"], row["thres"]))
}
  
create_predictions_boxplot <- function(dep_type, dep_level) {
  filename <- paste0("boxplots_pred_probs/", dep_type, "_", dep_level, "_preds_probs_boxplot.RDS")
  preds_tib <- weighted_preds_by_level(dep_type, dep_level) |> bind_rows() |>
    mutate(scenario = stringr::str_to_title(gsub("_", ", ", gsub(paste0(dep_type, "_", dep_level, "_"), "", scenario))))
  temp_plot <- preds_tib |> ggplot(aes(x = scenario, y = bma_preds, fill = bma_method)) + 
    geom_boxplot() +
    # geom_hline(yintercept = 8.210874e-05, col = "darkgrey", linetype = "longdash") + 
    theme_classic() +
    xlab("Likelihood and Threshold") + ylab("Prediction probabilities") + labs(fill = "")
  # ggsave(filename,
  #        plot = temp_plot,
  #        bg = 'transparent',
  #        dpi = 320)
  saveRDS(temp_plot, filename)
  print(paste0(filename, " has been saved"))
}

dep_types <- c("gauss", "logistic")
dep_levels <- c("high", "mid", "low")
types_levels_combos <- expand_grid(dep_types, dep_levels)
apply(types_levels_combos, 1, 
      function(row) create_predictions_boxplot(row["dep_types"], row["dep_levels"]))


# determining the true probability -----

lower_limit <- qnorm(pexp(8))
upper_limit <- qnorm(pexp(10))
pmvnorm(lower = rep(lower_limit, 2), upper = rep(upper_limit, 2), corr = matrix(c(1, 0.9, 0.9, 1), nrow = 2))[1]

lower_limit <- qgev(pexp(8), loc = 0, scale = 1, shape = 0)
upper_limit <- qgev(pexp(10), loc = 0, scale = 1, shape = 0)
pbvevd(q = c(lower_limit, upper_limit), dep = 0.1)

upper_rt_pt <- pbvevd(q = c(upper_limit, upper_limit), dep = 0.1)
upper_lt_pt <- pbvevd(q = c(lower_limit, upper_limit), dep = 0.1)
lower_rt_pt <- pbvevd(q = c(upper_limit, lower_limit), dep = 0.1)
lower_lt_pt <- pbvevd(q = c(lower_limit, lower_limit), dep = 0.1)

upper_rt_pt - upper_lt_pt - lower_rt_pt + lower_lt_pt
log_point <- rbvevd(5000, dep = 0.1)

## playing with using ALL angles -----------
med_pars <- gauss_high_trunc_marg_93[[1]][,1:2] |> as.numeric()
gw_fit <- gauss_gauge(W, 1-W, med_pars[2])
# gw_true <- gauss_gauge(W, 1-W, 0.9)
plot(W/gw_fit, (1-W)/gw_fit, col = 4)
# points(W/gw_true, (1-W)/gw_true, col = 3)

rw_gw_fit <- cbind(R, W, r0_w, gw_fit) |> as_tibble() |> mutate(r_gw = gw_fit * R) 
ctau <- quantile(rw_gw_fit$r_gw, 0.92)
rw_gw_fit <- rw_gw_fit |> mutate(r0w_tau = ctau / gw_fit)
rw_gw_fit_over1 <- rw_gw_fit |> filter(R > r0w_tau)

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
lower_limit <- qnorm(pexp(7))
upper_limit <- qnorm(pexp(9))
pmvnorm(lower = rep(lower_limit, 2), upper = rep(upper_limit, 2), corr = matrix(c(1, 0.5, 0.5, 1), nrow = 2))[1]

