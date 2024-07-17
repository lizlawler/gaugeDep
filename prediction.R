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

# importance weighting function
imp_weights <- function(k, w, r0w, pars, gauge) {
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  rate <- gauge_fcn(w, dep_par = pars[2:length(pars)])
  num <- pgamma(k * r0w, shape = pars[1], rate = rate, lower.tail = FALSE)
  denom <- pgamma(r0w, shape = pars[1], rate = rate, lower.tail = FALSE)
  return(num / denom)
}

# create function to simulate new data (bivariate case) -----------
sim_new_data <- function(k = 1, w, r0w, nsim, pars, gauge) {
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  if(k != 1) {
    weights <- imp_weights(k = k, w = w, r0w = r0w, pars = pars, gauge = gauge)
    resample_idx <- sample(1:length(w), size = nsim, prob = weights, replace = T)
    r0w_resample <- k * r0w[resample_idx]
  } else {
    resample_idx <- sample(1:length(w), size = nsim, replace = T)
    r0w_resample <- r0w[resample_idx]
  }
  w_resample <- w[resample_idx]
  rate_star <- gauge_fcn(w_resample, dep_par = pars[2:length(pars)])
  rstar <- qgamma(1 - runif(nsim) * pgamma(r0w_resample, shape = pars[1], 
                                           rate = rate_star, lower.tail = F),
                  shape = pars[1], rate = rate_star)
  xstar <- cbind(rstar * w_resample, rstar * (1 - w_resample))
  if(k != 1) {
    return(list(df = (xstar |> as_tibble() |> rename(x1 = V1, x2 = V2)),
                weights = weights))
  } else {
    return(list(df = (xstar |> as_tibble() |> rename(x1 = V1, x2 = V2))))
  }
}

# create function to make predictions from newly simulated data -------
pred_probs <- function(sim_df_list, idx, length_data, dim1, dim2, k = 1) {
  prob_x_given_r <- (sim_df_list[["df"]] |>
                       filter(x1 >= dim1[1] & x1 <= dim1[2] & x2 >= dim2[1] & x2 <= dim2[2]) |>
                       nrow()) / nrow(sim_df_list[["df"]])
  prob_r_over_thres <- length(idx)/length_data
  if(k != 1) {
    weights <- sim_df_list[["weights"]]
    prob_r_over_k_over_thres <- mean(weights)
    prob_r_over_k <- prob_r_over_k_over_thres * prob_r_over_thres
    return(prob_x_given_r * prob_r_over_k)
  } else {
    return(prob_x_given_r * prob_r_over_thres)
  }
}

# create function that wraps everything together ---------
make_preds <- function(data_file, posterior_pars, gauge, dim1, dim2, k = 1, all_angles = FALSE, true_gauge = FALSE) {
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  data <- fload(data_file)
  R <- data$R
  W <- data$W
  
  # use all angles or only angles associated with large R
  if(all_angles) {
    nsim <- 10000
    if(true_gauge) {
      gw_true <- data$ctau / data$r0_w_ctau
      ctau_low <- quantile(gw_true * R, 0.05)
      r0_w <- ctau_low / gw_true
    } else {
      gw_hat <- gauge_fcn(W, posterior_pars[2:length(posterior_pars)])
      ctau_hat <- quantile(R * gw_hat, 0.05)
      r0_w <- ctau_hat / gw_hat
    }
  } else {
    nsim <- 5000
    if(true_gauge) {
      r0_w <- data$r0_w_ctau
    } else {
      gw_hat <- gauge_fcn(W, posterior_pars[2:length(posterior_pars)])
      ctau_hat <- quantile(R * gw_hat, 0.95)
      r0_w <- ctau_hat / gw_hat
    }
  }
  idx <- which(R > r0_w)
  sim_df_list <- sim_new_data(k = k, w = W[idx], r0w = r0_w[idx], nsim = nsim, pars = posterior_pars, gauge = gauge)
  return(list(pred = pred_probs(sim_df_list, idx, length(W), dim1, dim2, k = k),
              new_data = sim_df_list))
}

# # testing -------
# gauss_params <- readRDS("extracted_params/gauss_gauss_high_cens_ctau_all_iter_params.RDS")
# rectangular_params <- readRDS("extracted_params/rectangular_gauss_low_cens_ctau_all_iter_params.RDS")
# median_params <- lapply(gauss_params, function(x) x |> select(-draw) |> apply(MARGIN = 2, FUN = median)) |> bind_rows()
# data_gauss_low_15 <- fload("data/gauss/low_15.json")
# b1_og <- make_preds("data/gauss/low_15.json",as.numeric(median_params[15,1:2]), gauge = "gauss", c(10,12), c(10,12), k = 1)
# b1 <- make_preds("data/gauss/low_15.json",as.numeric(median_params[15,1:2]), gauge = "gauss", c(10,12), c(10,12), k = 3.9)
# 
# b2 <- make_preds("data/gauss/low_15.json",as.numeric(median_params[15,1:2]), gauge = "gauss", c(10,12), c(6,8), k = 3.1)
# b3 <- make_preds("data/gauss/low_15.json",as.numeric(median_params[15,1:2]), gauge = "gauss", c(10,12), c(2,4), k = 2.4)
# 
# plot(b1$new_data$df, pch = 20)
# rect(xleft = 10, xright = 12, ybottom = 2, ytop = 4, border="red")
# rect(xleft = 10, xright = 12, ybottom = 6, ytop = 8, border="blue")
# rect(xleft = 10, xright = 12, ybottom = 10, ytop = 12, border="green")
# 
# gw_hat <- gauss_gauge(data_gauss_low_15$W, dep_par = 0.1)
# ctau_hat <- quantile(gw_hat * data_gauss_low_15$R, 0.95)
# test_grid <- expand_grid(x1_tilde = seq(10,12,length.out=25), x2_tilde = seq(6,8, length.out=25))
# test_grid <- test_grid |> mutate(w_tilde = x1_tilde / (x1_tilde + x2_tilde),
#                                  r_tilde = x1_tilde + x2_tilde)
# gw_tilde <- gauss_gauge(test_grid$w_tilde, 0.1)
# poss_k <- test_grid$r_tilde * gw_tilde / ctau_hat
# round(min(poss_k), 1) - 0.1


# create function to make predictions by the gauge function it was fit to ----------
preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, threshold, dim1, dim2, k = 1, all_angles = F, true_gauge = F) {
  posterior_params_all_iter <- readRDS(paste0("extracted_params/", gauge, "_", dep_type, "_", 
                                              dep_level, "_", likelihood, "_", threshold, "_all_iter_params.RDS"))
  init_data_path <- paste0("data/", dep_type, "/", dep_level, "_")
  posterior_params <- lapply(posterior_params_all_iter, function(x) x |> select(-draw) |> apply(MARGIN = 2, FUN = median)) |> bind_rows()
  results <- apply(posterior_params, 1, function(row) {
    data_file <- paste0(init_data_path, row["dataset"], ".json")
    params <- as.numeric(row[-length(row)])
    return(make_preds(data_file, posterior_pars = params, gauge = gauge, 
                      dim1 = dim1, dim2 = dim2, 
                      k = k, all_angles = all_angles, true_gauge = true_gauge))
  })
  return(lapply(results, function(x) x$pred) |> unlist())
}

# create function that makes predictions for all 100 datasets for a specific dependence type and level, likelihood type,
# threshold type, and with all gauge function fits -------
preds_by_dep_level_lhood_thres <- function(dep_type, dep_level, likelihood, threshold, dim1, dim2, k = 1, all_angles = F, true_gauge = F) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(sapply(gauge_library, function(x) preds_by_gauge(x, dep_type, dep_level, 
                                                          likelihood, threshold, 
                                                          dim1, dim2, 
                                                          k = k, all_angles = all_angles, true_gauge = true_gauge)) |>
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
weighted_preds_by_lhood_thres <- function(dep_type, dep_level, likelihood, threshold, dim1, dim2, k = 1, all_angles = F, true_gauge = F) {
  scenario <- paste0(dep_type, "_", dep_level, "_", likelihood, "_", threshold)
  temp_preds <- preds_by_dep_level_lhood_thres(dep_type, dep_level, 
                                               likelihood, threshold, 
                                               dim1, dim2, 
                                               k = k, all_angles = all_angles, true_gauge = true_gauge) |>
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

weighted_preds_by_level <- function(dep_type, dep_level, dim1, dim2, k = 1, all_angles = F, true_gauge = F) {
  thres <- c("ctau", "marg")
  lhood <- c("trunc", "cens")
  lhood_thres_combos <- expand_grid(lhood, thres)
  all_wts <- apply(lhood_thres_combos, 1, 
                   function(row) weighted_preds_by_lhood_thres(dep_type, dep_level, 
                                                               row["lhood"], row["thres"],
                                                               dim1, dim2, 
                                                               k = k, all_angles = all_angles, true_gauge = true_gauge))
}

# functions to determine the true probability -----
true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(pmvnorm(lower = c(dim1_star[1],dim2_star[1]), upper = c(dim1_star[2],dim2_star[2]), corr = corr_matrix)[1])
}

true_logistic_prob <- function(dim1, dim2, dep) {
  dim1_star <- qgev(pexp(dim1), loc = 0, scale = 1, shape = 0)
  dim2_star <- qgev(pexp(dim2), loc = 0, scale = 1, shape = 0)
  upper_right <- pbvevd(q = c(dim1_star[2], dim2_star[2]), dep = dep)
  upper_left <- pbvevd(q = c(dim1_star[1], dim2_star[2]), dep = dep)
  lower_right <- pbvevd(q = c(dim1_star[2], dim2_star[1]), dep = dep)
  lower_left <- pbvevd(q = c(dim1_star[1], dim2_star[1]), dep = dep)
  return(upper_right - upper_left - lower_right + lower_left)
}

# create function that wraps everything to gether to make boxplot of predictions
create_predictions_boxplot <- function(dep_type, dep_level, dim1, dim2, all_angles = F, true_gauge = F) {
  angles_name <- ifelse(all_angles, "all_angles", "angles_largeR")
  thresh_name <- ifelse(true_gauge, "true_gauge", "fitted_gauge")
  filename <- paste0("boxplots_pred_probs/", dep_type, "/", angles_name, "/",
                       thresh_name, "/", dep_level, "_preds_probs_boxplot.RDS")  
  
  # determine appropriate value of k
  pseudo_pred <- expand_grid(x1_pseudo = seq(dim1[1], dim1[2], length.out=25), 
                             x2_pseudo = seq(dim2[1], dim2[2], length.out=25)) |> 
    mutate(w_pseudo = x1_pseudo / (x1_pseudo + x2_pseudo),
           r_pseudo = x1_pseudo + x2_pseudo)
  if(dep_type == "gauss") {
    data_random <- RcppSimdJson::fload(paste0("data/gauss/", dep_level, "_", sample(1:100, size = 1),".json"))
    levels_list <- list(low = 0.1, mid = 0.5, high = 0.9, wc = 0.8)
    gw_pseudo <- gauss_gauge(pseudo_pred$w_pseudo, as.numeric(levels_list[dep_level]))
    if(all_angles) {
      ctau_temp <- data_random$ctau
      r0_w_temp <- data_random$r0_w_ctau
      gw_temp <- ctau_temp / r0_w_temp
      ctau_random <- quantile(gw_temp * data_random$R, 0.05)
    } else {
      ctau_random <- data$ctau
    }
    poss_k <- pseudo_pred$r_pseudo * gw_pseudo / ctau_random
    k <- max(round(min(poss_k), 1) - 0.2, 1)
  } else {
    data_random <- RcppSimdJson::fload(paste0("data/logistic/", dep_level, "_", sample(1:100, size = 1),".json"))
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1, wc_mid = 0.4, wc_low = 0.8)
    gw_pseudo <- logistic_gauge(pseudo_pred$w_pseudo, as.numeric(levels_list[dep_level]))
    if(all_angles) {
      ctau_temp <- data_random$ctau
      r0_w_temp <- data_random$r0_w_ctau
      gw_temp <- ctau_temp / r0_w_temp
      ctau_random <- quantile(gw_temp * data_random$R, 0.05)
    } else {
      ctau_random <- data$ctau
    }
    poss_k <- pseudo_pred$r_pseudo * gw_pseudo / ctau_random
    k <- max(round(min(poss_k), 1) - 0.1, 1)
  }
  
  plot_title <- paste0(dep_type, ", ", dep_level, ", (",paste(dim1, collapse = ","), ") x (", paste(dim2, collapse = ","),")", 
                       ", k = ", k, ", ", gsub("_", "", angles_name), ", threshold = ", gsub("_", " ", thresh_name))

  # determine true probability
  if(dep_type == "gauss") {
    levels_list <- list(low = 0.1, mid = 0.5, high = 0.9, wc = 0.8)
    true_prob <- true_gauss_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]))
  } else {
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1, wc_mid = 0.4, wc_low = 0.8)
    true_prob <- true_logistic_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]))
  }
  
  # make predictions
  preds_tib <- weighted_preds_by_level(dep_type, dep_level, 
                                       dim1, dim2, 
                                       k = k, all_angles = all_angles, true_gauge = true_gauge) |> 
    bind_rows() |>
    mutate(scenario = stringr::str_to_title(gsub("_", ", ", gsub(paste0(dep_type, "_", dep_level, "_"), "", scenario))))
  
  # create boxplot
  temp_plot <- preds_tib |> ggplot(aes(x = scenario, y = bma_preds, fill = bma_method)) + 
    geom_boxplot() +
    geom_hline(yintercept = true_prob, col = "darkgrey", linetype = "longdash") +
    theme_classic() +
    ggtitle(plot_title) +
    xlab("Likelihood and Threshold") + ylab("Prediction probabilities") + labs(fill = "")
  ggsave(gsub(".RDS", ".pdf", filename),
         plot = temp_plot,
         bg = 'transparent',
         dpi = 320)
  saveRDS(temp_plot, filename)
  print(paste0(filename, " has been saved"))
}

dep_types <- c("gauss")
dep_levels <- c("high", "mid", "low", "wc")
boxes <- tibble(dim1 = list(c(10,12)), dim2 = list(c(10,12), c(6,8), c(2,4)))
all_angles_vals <- c(T, F)
true_gauge_vals <- c(T, F)
all_combos_gauss <- expand_grid(dep_types, dep_levels, boxes, all_angles_vals, true_gauge_vals)
system.time(apply(all_combos_gauss, 1,
                  function(row) create_predictions_boxplot(dep_type = as.character(row["dep_types"]), 
                                                           dep_level = as.character(row["dep_levels"]),
                                                           dim1 = as.numeric(unlist(row["dim1"])), 
                                                           dim2 = as.numeric(unlist(row["dim2"])),
                                                           all_angles = as.logical(row["all_angles_vals"]), 
                                                           true_gauge = as.logical(row["true_gauge_vals"]))))


dep_types <- c("logistic")
dep_levels <- c("high", "mid", "low", "wc_mid", "wc_low")
all_combos_logistic <- expand_grid(dep_types, dep_levels, boxes, all_angles_vals, true_gauge_vals)
system.time(apply(all_combos_logistic, 1,
                  function(row) create_predictions_boxplot(dep_type = as.character(row["dep_types"]), 
                                                           dep_level = as.character(row["dep_levels"]),
                                                           dim1 = as.numeric(unlist(row["dim1"])), 
                                                           dim2 = as.numeric(unlist(row["dim2"])),
                                                           all_angles = as.logical(row["all_angles_vals"]), 
                                                           true_gauge = as.logical(row["true_gauge_vals"]))))
