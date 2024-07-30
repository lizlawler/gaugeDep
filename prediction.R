# library(geometricMVE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(progressr)
source("run_wc_models.R")

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

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

# w <- seq(0,1,length.out = 300)
# gw <- inv_log_gauge(w, dep_par = 31)
# plot(w/gw, (1-w)/gw, pch = 20)

# importance weighting function
imp_weights <- function(k, w, r0w, pars, gauge) {
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  rate <- gauge_fcn(w, dep_par = pars[2:length(pars)])
  num <- pgamma(k * r0w, shape = pars[1], rate = rate, lower.tail = FALSE, log.p = TRUE)
  denom <- pgamma(r0w, shape = pars[1], rate = rate, lower.tail = FALSE, log.p = TRUE)
  return(exp(num - denom))
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
make_preds <- function(data_file, posterior_pars, gauge, dim1, dim2, true_threshold = F, wc = F) {
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  data <- fload(data_file)
  R <- data$R
  W <- data$W
  ctau_true <- data$ctau
  r0_w_ctau <- data$r0_w_ctau
  
  # create fake data to use in determining k value
  pseudo_pred <- expand_grid(x1_pseudo = seq(dim1[1], dim1[2], length.out=15), 
                             x2_pseudo = seq(dim2[1], dim2[2], length.out=15)) |> 
    mutate(w_pseudo = x1_pseudo / (x1_pseudo + x2_pseudo),
           r_pseudo = x1_pseudo + x2_pseudo)
  
  if(true_threshold & !wc) {
    ro_w <- r0_w_ctau
  } else if(true_threshold & wc) {
    
  }
  # create threshold based on fitted gauge function (for LS and WC models)
  gw_fitted <- gauge_fcn(W, posterior_pars[2:length(posterior_pars)])
  ctau_fitted <- quantile(gw_fitted * R, 0.95)
  r0_w <- ctau_fitted/gw_fitted
  idx <- which(R > r0_w)
  
  # determine appropriate value of k with fitted gauge, in the proposed box
  gw_pseudo <- gauge_fcn(pseudo_pred$w_pseudo, posterior_pars[2:length(posterior_pars)])
  poss_k <- pseudo_pred$r_pseudo * gw_pseudo / ctau_fitted
  k <- max(round(min(poss_k), 1) - 0.1, 1)
  
  sim_df_list <- sim_new_data(k = k, w = W[idx], r0w = r0_w[idx], nsim = 5000, pars = posterior_pars, gauge = gauge)
  return(list(pred = pred_probs(sim_df_list, idx = idx, length_data = length(R), dim1, dim2, k = k),
              new_data = sim_df_list))
}

# create function to make predictions by the gauge function it was fit to ----------
preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, threshold, dim1, dim2) {
  posterior_params_all_iter <- readRDS(paste0("extracted_params/", gauge, "_", dep_type, "_", 
                                              dep_level, "_", likelihood, "_", threshold, "_all_iter_params.RDS"))
  init_data_path <- paste0("data/", dep_type, "/", dep_level, "_")
  posterior_params <- lapply(posterior_params_all_iter, function(x) x |> select(-draw) |> apply(MARGIN = 2, FUN = median)) |> bind_rows()
  results <- apply(posterior_params, 1, function(row) {
    data_file <- paste0(init_data_path, row["dataset"], ".json")
    params <- as.numeric(row[-length(row)])
    return(tryCatch(make_preds(data_file, posterior_pars = params, gauge = gauge, 
                      dim1 = dim1, dim2 = dim2), error=function(e) list(pred = NA)))
  })
  return(lapply(results, function(x) x$pred) |> unlist())
}

# create function that makes predictions for all 100 datasets for a specific dependence type and level, likelihood type,
# threshold type, and with all gauge function fits -------
preds_by_dep_level_lhood_thres <- function(dep_type, dep_level, likelihood, threshold, dim1, dim2) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(sapply(gauge_library, function(x) preds_by_gauge(x, dep_type, dep_level, 
                                                          likelihood, threshold, 
                                                          dim1, dim2)) |>
           as_tibble() |>
           mutate(dataset = 1:100))
}

# test <- preds_by_dep_level_lhood_thres("gauss", "high", "trunc", "ctau", c(10,12), c(10,12))
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
weighted_preds_by_lhood_thres <- function(dep_type, dep_level, likelihood, threshold, dim1, dim2) {
  scenario <- paste0(dep_type, "_", dep_level, "_", likelihood, "_", threshold)
  temp_preds <- preds_by_dep_level_lhood_thres(dep_type, dep_level, 
                                               likelihood, threshold, 
                                               dim1, dim2) |>
    pivot_longer(cols = -'dataset', names_to = "method", values_to = "preds")
  temp_wts <- make_wts_df(paste0("stacking_weights/", scenario, "_wts.RDS"))
  temp_weighted_preds <- suppressMessages(temp_wts |> left_join(temp_preds) |>
                                            mutate(stacking_preds = preds * stacking,
                                                   pseudo_boot = pseudobma_boot * preds,
                                                   pseudo_noboot = pseudobma_noboot * preds) |>
                                            group_by(dataset) |>
                                            summarize(stacking_predictions = sum(stacking_preds),
                                                      pseudobma_boot_preds = sum(pseudo_boot),
                                                      pseudobma_noboot_preds = sum(pseudo_noboot)) |>
                                            ungroup())
  boxplot_wts <- temp_weighted_preds |> 
    pivot_longer(cols = -'dataset', names_to = "method", values_to = "preds") |>
    mutate(method = case_when(grepl("stacking", method) ~ 'Stacking',
                              grepl("noboot", method) ~ 'Pseudo-BMA',
                              grepl("boot", method) ~ 'Pseudo-BMA+'),
           method = as.factor(method),
           scenario = rep(scenario))
  return(boxplot_wts)
}

weighted_preds_by_level <- function(dep_type, dep_level, dim1, dim2) {
  thres <- c("ctau", "marg")
  lhood <- c("trunc", "cens")
  lhood_thres_combos <- expand_grid(lhood, thres)
  all_wts <- apply(lhood_thres_combos, 1, function(row) {
    weighted_preds_by_lhood_thres(dep_type, dep_level, 
                                  row["lhood"], row["thres"],
                                  dim1, dim2)})
  return(all_wts)
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
create_predictions_boxplot <- function(dep_type, dep_level, box_num) {
  dim1 <- c(10, 12)
  if(box_num == "b1") {
    dim2 <- dim1
  } else if(box_num == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  plot_filename <- paste0("boxplots_pred_probs/", dep_type, "/", 
                          dep_level, "_", box_num, "_preds_boxplot_with_wc.pdf") 
  rds_filename <- paste0("boxplots_pred_probs/", dep_type, "/rds_files/", dep_level, 
                         "_",  box_num, "_preds_boxplot_with_wc.RDS")
  
  # make predictions for Lawer and Shaby method
  ls_preds_tib <- weighted_preds_by_level(dep_type, dep_level, 
                                          dim1, dim2) |> 
    bind_rows() |>
    mutate(scenario = stringr::str_to_title(gsub("_", ", ", gsub(paste0(dep_type, "_", dep_level, "_"), "", scenario))))
  # 
  # make predictions for Wadsworth and Campbell method
  wc_preds_tib <- apply(split_wc_fits[[paste0(dep_type,".",dep_level)]], 1,
                        function(row) tryCatch(make_preds(paste0("data/", dep_type, "/", dep_level, "_", row["datasets"], ".json"),
                                     posterior_pars = as.numeric(unlist(row["mle"])),
                                     gauge = row["gauge_name"],
                                     dim1 = dim1,
                                     dim2 = dim2), error=function(e) list(pred=NA))) |>
                            lapply(function(x) x$pred) |>
                            unlist() |>
                            as_tibble() |>
                            rename(preds = value) |>
                            mutate(dataset = 1:100, scenario = 'W-C', method = NA)

  all_preds <- ls_preds_tib |> rbind(wc_preds_tib)
  
  # determine true probability
  if(dep_type == "gauss") {
    levels_list <- list(low = 0.1, mid = 0.5, high = 0.9)
    true_prob <- true_gauss_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]))
  } else {
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1)
    true_prob <- true_logistic_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]))
  }
  
  # create boxplot
  plot_title <- paste0(dep_type, ", ", dep_level, ", (",paste(dim1, collapse = ","), ") x (", paste(dim2, collapse = ","),")")
  temp_plot <- all_preds |> ggplot(aes(x = scenario, y = preds, fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = true_prob, col = "darkgrey", linetype = "longdash") +
    theme_classic() +
    ggtitle(plot_title) +
    xlab("Likelihood and Threshold") + ylab("Prediction probabilities") + labs(fill = "")
  ggsave(plot_filename,
         plot = temp_plot,
         bg = 'transparent',
         width = 8,
         height = 7,
         dpi = 320)
  saveRDS(temp_plot, rds_filename)
  print(paste0(plot_filename, " has been saved"))
}

dep_types <- c("gauss", "logistic")
dep_levels <- c("high", "mid", "low")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)
all_combos <- all_combos[12:18,]

with_progress({
  # Create a progress handler
  p <- progressor(steps = nrow(all_combos))
  
  # Apply the function using apply and update the progress bar
  apply(all_combos, 1, function(row) {
    p()  # Update the progress bar
    create_predictions_boxplot(dep_type = row["dep_types"], 
                               dep_level = row["dep_levels"],
                               box_num = row["boxes"])
  })
})

