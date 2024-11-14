# library(geometricMVE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(progressr)
source("gauge_functions_wrt_w.R")
source("run_wc_models.R")

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

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
make_preds <- function(data_file, posterior_pars, gauge, dim1, dim2, true_threshold = F, wc = F, k_var = T) {
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  data <- fload(data_file)
  R <- data$R
  W <- data$W
  
  # determine threshold for the fitted gauge with posterior (or MLE) parameters
  gw_fitted <- gauge_fcn(W, posterior_pars[2:length(posterior_pars)])
  ctau_fitted <- quantile(gw_fitted * R, 0.95)
  
  if(k_var) {
    # create fake data to use in determining k value
    pseudo_pred <- expand_grid(x1_pseudo = seq(dim1[1], dim1[2], length.out=15), 
                               x2_pseudo = seq(dim2[1], dim2[2], length.out=15)) |> 
      mutate(w_pseudo = x1_pseudo / (x1_pseudo + x2_pseudo),
             r_pseudo = x1_pseudo + x2_pseudo)
    # determine ideal value of k using the above
    gw_pseudo <- gauge_fcn(pseudo_pred$w_pseudo, posterior_pars[2:length(posterior_pars)])
    poss_k <- pseudo_pred$r_pseudo * gw_pseudo / ctau_fitted
    k <- max(round(min(poss_k), 1) - 0.1, 1)
  } else {
    k <- 1
  }
  
  if(true_threshold & !wc) {
    # pull true gauge function threshold from data
    r0_w <- data$r0_w_ctau
  } else if (true_threshold & wc) {
    # recreate empirical threshold
    temp_qr <- geometricMVE::QR.2d(r = R, w = W, method = "empirical")
    r0_w <- temp_qr$r0w     
  } else {
    # create threshold based on fitted gauge function (for LS and WC models)
    r0_w <- ctau_fitted/gw_fitted
  }
  
  idx <- which(R > r0_w)
  sim_df_list <- sim_new_data(k = k, w = W[idx], r0w = r0_w[idx], nsim = 5000, pars = posterior_pars, gauge = gauge)
  return(list(pred = pred_probs(sim_df_list, idx = idx, length_data = length(R), dim1, dim2, k = k),
              new_data = sim_df_list))
}

# create function to make predictions by the gauge function it was fit to ----------
preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, threshold, dim1, dim2, 
                           true_threshold = F, wc = F, k_var = T) {
  posterior_params_all_iter <- readRDS(paste0("extracted_params/", gauge, "_", dep_type, "_", 
                                              dep_level, "_", likelihood, "_", threshold, "_all_iter_params.RDS"))
  init_data_path <- paste0("data/", dep_type, "/", dep_level, "_")
  posterior_params <- lapply(posterior_params_all_iter, function(x) x |> select(-draw) |> apply(MARGIN = 2, FUN = median)) |> bind_rows()
  results <- apply(posterior_params, 1, function(row) {
    data_file <- paste0(init_data_path, row["dataset"], ".json")
    params <- as.numeric(row[-length(row)])
    return(tryCatch(make_preds(data_file, posterior_pars = params, gauge = gauge, 
                               dim1 = dim1, dim2 = dim2,
                               true_threshold = true_threshold, wc = wc, k_var = k_var), 
                    error=function(e) list(pred = NA)))
  })
  return(lapply(results, function(x) x$pred) |> unlist())
}

# create function that makes predictions for all 100 datasets for a specific dependence type and level, likelihood type,
# threshold type, and with all gauge function fits -------
preds_by_dep_level_lhood_thres <- function(dep_type, dep_level, likelihood, threshold, dim1, dim2, 
                                           true_threshold = F, wc = F, k_var = T) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(sapply(gauge_library, function(x) preds_by_gauge(x, dep_type, dep_level, 
                                                          likelihood, threshold, 
                                                          dim1, dim2,
                                                          true_threshold, wc, k_var)) |>
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
weighted_preds_by_lhood_thres <- function(dep_type, dep_level, likelihood, threshold, dim1, dim2, 
                                          true_threshold = F, wc = F, k_var = T) {
  scenario <- paste0(dep_type, "_", dep_level, "_", likelihood, "_", threshold)
  temp_preds <- preds_by_dep_level_lhood_thres(dep_type, dep_level, 
                                               likelihood, threshold, 
                                               dim1, dim2,
                                               true_threshold, wc, k_var) |>
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

weighted_preds_by_level <- function(dep_type, dep_level, dim1, dim2, true_threshold = F, wc = F, k_var = T) {
  thres <- c("ctau", "marg")
  lhood <- c("trunc", "cens")
  lhood_thres_combos <- expand_grid(lhood, thres)
  if(dep_type == "husler_reiss") {
    lhood_thres_combos <- lhood_thres_combos |> filter(thres != "ctau")
  }
  all_wts <- apply(lhood_thres_combos, 1, function(row) {
    weighted_preds_by_lhood_thres(dep_type, dep_level, 
                                  row["lhood"], row["thres"],
                                  dim1, dim2,
                                  true_threshold, wc, k_var)})
  return(all_wts)
}

# functions to determine the true probability -----
true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(pmvnorm(lower = c(dim1_star[1],dim2_star[1]), upper = c(dim1_star[2],dim2_star[2]), corr = corr_matrix)[1])
}

true_bvevd_prob <- function(dim1, dim2, dep, model_type) {
  dim1_star <- qgev(pexp(dim1), loc = 0, scale = 1, shape = 0)
  dim2_star <- qgev(pexp(dim2), loc = 0, scale = 1, shape = 0)
  upper_right <- pbvevd(q = c(dim1_star[2], dim2_star[2]), model = model_type, dep = dep)
  upper_left <- pbvevd(q = c(dim1_star[1], dim2_star[2]), model = model_type, dep = dep)
  lower_right <- pbvevd(q = c(dim1_star[2], dim2_star[1]), model = model_type, dep = dep)
  lower_left <- pbvevd(q = c(dim1_star[1], dim2_star[1]), model = model_type, dep = dep)
  return(upper_right - upper_left - lower_right + lower_left)
}

# create function that wraps everything to gether to make boxplot of predictions
create_predictions_boxplot <- function(dep_type, dep_level, box_num, true_threshold = F, wc = F, k_var = T) {
  dim1 <- c(10, 12)
  if(box_num == "b1") {
    dim2 <- dim1
  } else if(box_num == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  thresh_name <- ifelse(true_threshold, "true", "fitted")
  k_file <- ifelse(k_var, "k_not1", "k_is1")
  k_title <- ifelse(k_var, "k > 1", "k = 1")
  plot_filename <- paste0("boxplots_pred_probs/", dep_type, "/", 
                          dep_level, "_", box_num, "_", k_file, "_",
                          thresh_name, "_threshold_preds_boxplot_with_wc.pdf") 
  rds_filename <- paste0("boxplots_pred_probs/", dep_type, "/rds_files/", 
                         dep_level, "_", box_num, "_", k_file, "_",
                         thresh_name, "_threshold_preds_boxplot_with_wc.RDS")
  
  # make predictions for Lawer and Shaby method
  ls_preds_tib <- weighted_preds_by_level(dep_type, dep_level, 
                                          dim1, dim2,
                                          true_threshold, wc = F, k_var) |> 
    bind_rows() |>
    mutate(scenario = stringr::str_to_title(gsub("_", ", ", gsub(paste0(dep_type, "_", dep_level, "_"), "", scenario))))
  # 
  # make predictions for Wadsworth and Campbell method
  wc_preds_tib <- apply(split_wc_fits[[paste0(dep_type,".",dep_level)]], 1,
                        function(row) tryCatch(make_preds(paste0("data/", dep_type, "/", dep_level, "_", row["datasets"], ".json"),
                                                          posterior_pars = as.numeric(unlist(row["mle"])),
                                                          gauge = row["gauge_name"],
                                                          dim1 = dim1,
                                                          dim2 = dim2,
                                                          true_threshold, wc = T, k_var), 
                                               error=function(e) list(pred=NA))) |>
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
  } else if(dep_type == "logistic"){
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "log")
  } else {
    levels_list <- list(low = 0.25, mid = 2, high = 6)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "hr")
  }
  # create boxplot
  plot_title <- paste0(dep_type, ", ", dep_level, 
                       ", (",paste(dim1, collapse = ","), ") x (", paste(dim2, collapse = ","),"), ", 
                       k_title, ", ", thresh_name, " threshold")
  temp_plot <- all_preds |> ggplot(aes(x = scenario, y = preds, fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = true_prob, col = "darkgrey", linetype = "longdash") +
    theme_classic() +
    ggtitle(plot_title) + scale_fill_discrete(breaks = ~ .x[!is.na(.x)]) +
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

# create_predictions_boxplot("gauss", "high", "b3", T, F)
dep_types <- c("gauss", "logistic", "husler_reiss")
dep_levels <- c("high", "mid", "low")
boxes <- c("b1", "b2", "b3")
true_thresh_vals <- c(T, F)
k_vals <- c(T, F)
all_combos <- expand_grid(dep_types, dep_levels, boxes, true_thresh_vals, k_vals) |>
  filter(!(dep_types == "husler_reiss" & true_thresh_vals == T)) # no "true gauge" threshold for HR case

with_progress({
  # Create a progress handler
  p <- progressor(steps = nrow(all_combos))
  
  # Apply the function using apply and update the progress bar
  apply(all_combos, 1, function(row) {
    p()  # Update the progress bar
    create_predictions_boxplot(dep_type = row["dep_types"], 
                               dep_level = row["dep_levels"],
                               box_num = row["boxes"], 
                               true_threshold = as.logical(row["true_thresh_vals"]),
                               wc = F,
                               k_var = as.logical(row["k_vals"]))
  })
})


# function to read each plot file and appropriately rename
read_files <- function(file, plot_name) {
  temp <- readRDS(file)
  assign(plot_name, temp, parent.frame())
  rm(temp)
  gc()
}
## read in Gaussian boxplots and patch together -------
all_gauss_files <- paste0("boxplots_pred_probs/gauss/rds_files/",
                          list.files("boxplots_pred_probs/gauss/rds_files/",
                                     pattern = "threshold_preds"))
gauss_true_files <- all_gauss_files[grepl("true", all_gauss_files)]
gauss_fitted_files <- all_gauss_files[grepl("fitted", all_gauss_files)]
gauss_true_k_not1_files <- gauss_true_files[grepl("not1", gauss_true_files)]
gauss_true_k_not1_names <- str_remove(basename(gauss_true_k_not1_files), "_threshold_preds_boxplot_with_wc.RDS")
gauss_true_k_is1_files <- gauss_true_files[grepl("is1", gauss_true_files)]
gauss_true_k_is1_names <- str_remove(basename(gauss_true_k_is1_files), "_threshold_preds_boxplot_with_wc.RDS")

gauss_fitted_k_not1_files <- gauss_fitted_files[grepl("not1", gauss_fitted_files)]
gauss_fitted_k_not1_names <- str_remove(basename(gauss_fitted_k_not1_files), "_threshold_preds_boxplot_with_wc.RDS")
gauss_fitted_k_is1_files <- gauss_fitted_files[grepl("is1", gauss_fitted_files)]
gauss_fitted_k_is1_names <- str_remove(basename(gauss_fitted_k_is1_files), "_threshold_preds_boxplot_with_wc.RDS")


for(i in seq_along(gauss_true_k_not1_files)) {
  read_files(gauss_true_k_not1_files[i], gauss_true_k_not1_names[i])
}

# for(i in seq_along(gauss_true_k_is1_files)) {
#   read_files(gauss_true_k_is1_files[i], gauss_true_k_is1_names[i])
# }

for(i in seq_along(gauss_fitted_k_not1_files)) {
  read_files(gauss_fitted_k_not1_files[i], gauss_fitted_k_not1_names[i])
}

# for(i in seq_along(gauss_fitted_k_is1_files)) {
#   read_files(gauss_fitted_k_is1_files[i], gauss_fitted_k_is1_names[i])
# }

gauss_high_b1_axis <- range(ggplot_build(high_b1_k_not1_true)$layout$panel_params[[1]]$y.range)
gauss_high_b2_axis <- range(ggplot_build(high_b2_k_not1_true)$layout$panel_params[[1]]$y.range)
gauss_high_b3_axis <- range(ggplot_build(high_b3_k_not1_true)$layout$panel_params[[1]]$y.range)

gauss_mid_b1_axis <- range(ggplot_build(mid_b1_k_not1_true)$layout$panel_params[[1]]$y.range)
gauss_mid_b2_axis <- range(ggplot_build(mid_b2_k_not1_true)$layout$panel_params[[1]]$y.range)
gauss_mid_b3_axis <- range(ggplot_build(mid_b3_k_not1_true)$layout$panel_params[[1]]$y.range)

gauss_low_b1_axis <- range(ggplot_build(low_b1_k_not1_true)$layout$panel_params[[1]]$y.range)
gauss_low_b2_axis <- range(ggplot_build(low_b2_k_not1_true)$layout$panel_params[[1]]$y.range)
gauss_low_b3_axis <- range(ggplot_build(low_b3_k_not1_true)$layout$panel_params[[1]]$y.range)

## create plots with all 3 boxes together (Gauss) ---------------
# using true threshold (Gauss) ------
library(patchwork)
gauss_high_true_k_not1_all <- (high_b1_k_not1_true + ggtitle(NULL)) + 
  (high_b2_k_not1_true + ggtitle(NULL)) + 
  (high_b3_k_not1_true + ggtitle(NULL) + coord_cartesian(ylim = c(0, 2.5e-7))) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
gauss_high_b3_axis <- c(0, 2.5e-7)

ggsave("bma_update_deck/gauss_high_true_k_not1.pdf",
       plot = gauss_high_true_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# gauss_high_true_k_is1_all <- (high_b1_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b1_axis)) + 
#   (high_b2_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b2_axis)) + 
#   (high_b3_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = c(0, 2.5e-7))) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))

gauss_mid_true_k_not1_all <- (mid_b1_k_not1_true + ggtitle(NULL)) + 
  (mid_b2_k_not1_true + ggtitle(NULL)) + 
  (mid_b3_k_not1_true + ggtitle(NULL)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("bma_update_deck/gauss_mid_true_k_not1.pdf",
       plot = gauss_mid_true_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# gauss_mid_true_k_is1_all <- (mid_b1_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b1_axis)) + 
#   (mid_b2_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b2_axis)) + 
#   (mid_b3_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# 
# ggsave("bma_update_deck/gauss_mid_true_k_is1.pdf",
#        plot = gauss_mid_true_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

gauss_low_true_k_not1_all <- (low_b1_k_not1_true + ggtitle(NULL) + coord_cartesian(ylim = c(0, 1e-7))) + 
  (low_b2_k_not1_true + ggtitle(NULL)) + 
  (low_b3_k_not1_true + ggtitle(NULL)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
gauss_low_b1_axis <- c(0, 1e-7)
ggsave("bma_update_deck/gauss_low_true_k_not1.pdf",
       plot = gauss_low_true_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# gauss_low_true_k_is1_all <- (low_b1_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b1_axis)) + 
#   (low_b2_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b2_axis)) + 
#   (low_b3_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'), 
#         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))

# ggsave("bma_update_deck/gauss_low_true_k_is1.pdf",
#        plot = gauss_low_true_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

## using fitted threshold (Gauss) --------
gauss_low_fitted_k_not1_all <- (low_b1_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b1_axis)) + 
  (low_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b2_axis)) + 
  (low_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/gauss_low_fitted_k_not1.pdf",
       plot = gauss_low_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# gauss_low_fitted_k_is1_all <- (low_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b1_axis)) + 
#   (low_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b2_axis)) + 
#   (low_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_low_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/gauss_low_fitted_k_is1.pdf",
#        plot = gauss_low_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

gauss_mid_fitted_k_not1_all <- (mid_b1_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b1_axis)) + 
  (mid_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b2_axis)) + 
  (mid_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/gauss_mid_fitted_k_not1.pdf",
       plot = gauss_mid_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# gauss_mid_fitted_k_is1_all <- (mid_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b1_axis)) + 
#   (mid_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b2_axis)) + 
#   (mid_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_mid_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/gauss_mid_fitted_k_is1.pdf",
#        plot = gauss_mid_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

gauss_high_fitted_k_not1_all <- (high_b1_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b1_axis)) + 
  (high_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b2_axis)) + 
  (high_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/gauss_high_fitted_k_not1.pdf",
       plot = gauss_high_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# gauss_high_fitted_k_is1_all <- (high_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b1_axis)) + 
#   (high_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b2_axis)) + 
#   (high_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = gauss_high_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/gauss_high_fitted_k_is1.pdf",
#        plot = gauss_high_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

## read in logistic boxplots and patch together -------
all_logistic_files <- paste0("boxplots_pred_probs/logistic/rds_files/",
                          list.files("boxplots_pred_probs/logistic/rds_files/",
                                     pattern = "threshold_preds"))
logistic_true_files <- all_logistic_files[grepl("true", all_logistic_files)]
logistic_fitted_files <- all_logistic_files[grepl("fitted", all_logistic_files)]
logistic_true_k_not1_files <- logistic_true_files[grepl("not1", logistic_true_files)]
logistic_true_k_not1_names <- str_remove(basename(logistic_true_k_not1_files), "_threshold_preds_boxplot_with_wc.RDS")
# logistic_true_k_is1_files <- logistic_true_files[grepl("is1", logistic_true_files)]
# logistic_true_k_is1_names <- str_remove(basename(logistic_true_k_is1_files), "_threshold_preds_boxplot_with_wc.RDS")

logistic_fitted_k_not1_files <- logistic_fitted_files[grepl("not1", logistic_fitted_files)]
logistic_fitted_k_not1_names <- str_remove(basename(logistic_fitted_k_not1_files), "_threshold_preds_boxplot_with_wc.RDS")
logistic_fitted_k_is1_files <- logistic_fitted_files[grepl("is1", logistic_fitted_files)]
logistic_fitted_k_is1_names <- str_remove(basename(logistic_fitted_k_is1_files), "_threshold_preds_boxplot_with_wc.RDS")

for(i in seq_along(logistic_true_k_not1_files)) {
  read_files(logistic_true_k_not1_files[i], logistic_true_k_not1_names[i])
}

# for(i in seq_along(logistic_true_k_is1_files)) {
#   read_files(logistic_true_k_is1_files[i], logistic_true_k_is1_names[i])
# }

for(i in seq_along(logistic_fitted_k_not1_files)) {
  read_files(logistic_fitted_k_not1_files[i], logistic_fitted_k_not1_names[i])
}

# for(i in seq_along(logistic_fitted_k_is1_files)) {
#   read_files(logistic_fitted_k_is1_files[i], logistic_fitted_k_is1_names[i])
# }

# pull axis limits so various scenarios are comparable
logistic_high_b1_axis <- range(ggplot_build(high_b1_k_not1_true)$layout$panel_params[[1]]$y.range)
logistic_high_b2_axis <- range(ggplot_build(high_b2_k_not1_true)$layout$panel_params[[1]]$y.range)
logistic_high_b3_axis <- range(ggplot_build(high_b3_k_not1_true)$layout$panel_params[[1]]$y.range)

logistic_mid_b1_axis <- range(ggplot_build(mid_b1_k_not1_true)$layout$panel_params[[1]]$y.range)
logistic_mid_b2_axis <- range(ggplot_build(mid_b2_k_not1_true)$layout$panel_params[[1]]$y.range)
logistic_mid_b3_axis <- range(ggplot_build(mid_b3_k_not1_true)$layout$panel_params[[1]]$y.range)

logistic_low_b1_axis <- range(ggplot_build(low_b1_k_not1_true)$layout$panel_params[[1]]$y.range)
logistic_low_b2_axis <- range(ggplot_build(low_b2_k_not1_true)$layout$panel_params[[1]]$y.range)
logistic_low_b3_axis <- range(ggplot_build(low_b3_k_not1_true)$layout$panel_params[[1]]$y.range)

## create plots with all 3 boxes together (logistic) ---------------
# using true threshold (logistic) -------
logistic_high_true_k_not1_all <- (high_b1_k_not1_true + ggtitle(NULL)) + 
  (high_b2_k_not1_true + ggtitle(NULL) + coord_cartesian(ylim = c(0, 5e-12))) + 
  (high_b3_k_not1_true + ggtitle(NULL) + coord_cartesian(ylim = c(0, 5e-14))) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

logistic_high_b3_axis <- c(0, 5e-14)
logistic_high_b2_axis <- c(0, 5e-12)

ggsave("bma_update_deck/logistic_high_true_k_not1.pdf",
       plot = logistic_high_true_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# logistic_high_true_k_is1_all <- (high_b1_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b1_axis)) + 
#   (high_b2_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b2_axis)) + 
#   (high_b3_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# 
# ggsave("bma_update_deck/logistic_high_true_k_is1_all.pdf",
#        plot = logistic_high_true_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

logistic_mid_true_k_not1_all <- (mid_b1_k_not1_true + ggtitle(NULL)) + 
  (mid_b2_k_not1_true + ggtitle(NULL)) + 
  (mid_b3_k_not1_true + ggtitle(NULL)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("bma_update_deck/logistic_mid_true_k_not1.pdf",
       plot = logistic_mid_true_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# logistic_mid_true_k_is1_all <- (mid_b1_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b1_axis)) + 
#   (mid_b2_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b2_axis)) + 
#   (mid_b3_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# 
# ggsave("bma_update_deck/logistic_mid_true_k_is1.pdf",
#        plot = logistic_mid_true_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

logistic_low_true_k_not1_all <- (low_b1_k_not1_true + ggtitle(NULL)) + 
  (low_b2_k_not1_true + ggtitle(NULL)) + 
  (low_b3_k_not1_true + ggtitle(NULL)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/logistic_low_true_k_not1.pdf",
       plot = logistic_low_true_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# logistic_low_true_k_is1_all <- (low_b1_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b1_axis)) + 
#   (low_b2_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b2_axis)) + 
#   (low_b3_k_is1_true + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# 
# ggsave("bma_update_deck/logistic_low_true_k_is1.pdf",
#        plot = logistic_low_true_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

## using fitted threshold (logistic) ---------
logistic_low_fitted_k_not1_all <- (low_b1_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b1_axis)) + 
  (low_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b2_axis)) + 
  (low_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/logistic_low_fitted_k_not1.pdf",
       plot = logistic_low_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# logistic_low_fitted_k_is1_all <- (low_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b1_axis)) + 
#   (low_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b2_axis)) + 
#   (low_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_low_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/logistic_low_fitted_k_is1.pdf",
#        plot = logistic_low_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

logistic_mid_fitted_k_not1_all <- (mid_b1_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b1_axis)) + 
  (mid_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b2_axis)) + 
  (mid_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/logistic_mid_fitted_k_not1.pdf",
       plot = logistic_mid_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# logistic_mid_fitted_k_is1_all <- (mid_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b1_axis)) + 
#   (mid_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b2_axis)) + 
#   (mid_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_mid_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         
#         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/logistic_mid_fitted_k_is1.pdf",
#        plot = logistic_mid_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

logistic_high_fitted_k_not1_all <- (high_b1_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b1_axis)) + 
  (high_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b2_axis)) + 
  (high_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/logistic_high_fitted_k_not1.pdf",
       plot = logistic_high_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# logistic_high_fitted_k_is1_all <- (high_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b1_axis)) + 
#   (high_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b2_axis)) + 
#   (high_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = logistic_high_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/logistic_high_fitted_k_is1.pdf",
#        plot = logistic_high_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)


## read in Husler Reiss boxplots and patch together ------------
all_husler_reiss_files <- paste0("boxplots_pred_probs/husler_reiss/rds_files/",
                             list.files("boxplots_pred_probs/husler_reiss/rds_files/",
                                        pattern = "threshold_preds"))
husler_reiss_fitted_k_not1_files <- all_husler_reiss_files[grepl("not1", all_husler_reiss_files)]
husler_reiss_fitted_k_not1_names <- str_remove(basename(husler_reiss_fitted_k_not1_files), "_threshold_preds_boxplot_with_wc.RDS")
husler_reiss_fitted_k_is1_files <- all_husler_reiss_files[grepl("is1", all_husler_reiss_files)]
husler_reiss_fitted_k_is1_names <- str_remove(basename(husler_reiss_fitted_k_is1_files), "_threshold_preds_boxplot_with_wc.RDS")

for(i in seq_along(husler_reiss_fitted_k_not1_files)) {
  read_files(husler_reiss_fitted_k_not1_files[i], husler_reiss_fitted_k_not1_names[i])
}

for(i in seq_along(husler_reiss_fitted_k_is1_files)) {
  read_files(husler_reiss_fitted_k_is1_files[i], husler_reiss_fitted_k_is1_names[i])
}

# pull axis limits so various scenarios are comparable
husler_reiss_high_b1_axis <- range(ggplot_build(high_b1_k_not1_fitted)$layout$panel_params[[1]]$y.range)
husler_reiss_high_b2_axis <- range(ggplot_build(high_b2_k_not1_fitted)$layout$panel_params[[1]]$y.range)
husler_reiss_high_b3_axis <- range(ggplot_build(high_b3_k_not1_fitted)$layout$panel_params[[1]]$y.range)

husler_reiss_mid_b1_axis <- range(ggplot_build(mid_b1_k_not1_fitted)$layout$panel_params[[1]]$y.range)
husler_reiss_mid_b2_axis <- range(ggplot_build(mid_b2_k_not1_fitted)$layout$panel_params[[1]]$y.range)
husler_reiss_mid_b3_axis <- range(ggplot_build(mid_b3_k_not1_fitted)$layout$panel_params[[1]]$y.range)

husler_reiss_low_b1_axis <- range(ggplot_build(low_b1_k_not1_fitted)$layout$panel_params[[1]]$y.range)
husler_reiss_low_b2_axis <- range(ggplot_build(low_b2_k_not1_fitted)$layout$panel_params[[1]]$y.range)
husler_reiss_low_b3_axis <- range(ggplot_build(low_b3_k_not1_fitted)$layout$panel_params[[1]]$y.range)

## create plots with all 3 boxes together (husler_reiss) ---------------
# using true threshold (husler_reiss) -------
husler_reiss_high_fitted_k_not1_all <- (high_b1_k_not1_fitted + ggtitle(NULL)) + 
  (high_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = c(0, 2.5e-9))) + 
  (high_b3_k_not1_fitted + ggtitle(NULL)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
husler_reiss_high_b2_axis <- c(0, 2.5e-9)
ggsave("bma_update_deck/husler_reiss_high_fitted_k_not1.pdf",
       plot = husler_reiss_high_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# husler_reiss_high_fitted_k_is1_all <- (high_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_high_b1_axis)) + 
#   (high_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_high_b2_axis)) + 
#   (high_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_high_b3_axis)) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/husler_reiss_high_fitted_k_is1_all.pdf",
#        plot = husler_reiss_high_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

husler_reiss_mid_b3_axis <- c(0, 5e-8)
husler_reiss_mid_fitted_k_not1_all <- (mid_b1_k_not1_fitted + ggtitle(NULL)) + 
  (mid_b2_k_not1_fitted + ggtitle(NULL)) + 
  (mid_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_mid_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/husler_reiss_mid_fitted_k_not1.pdf",
       plot = husler_reiss_mid_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

# husler_reiss_mid_fitted_k_is1_all <- (mid_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_mid_b1_axis)) + 
#   (mid_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_mid_b2_axis)) + 
#   (mid_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = c(0, 5e-12))) +
#   plot_layout(guides = 'collect') & 
#   theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# ggsave("bma_update_deck/husler_reiss_mid_fitted_k_is1_all.pdf",
#        plot = husler_reiss_mid_fitted_k_is1_all,
#        dpi = 320,
#        bg = "transparent",
#        width = 16, height =7)

husler_reiss_low_fitted_k_not1_all <- (low_b1_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = c(0, 1e-7))) + 
  (low_b2_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = c(0, 2.5e-6))) + 
  (low_b3_k_not1_fitted + ggtitle(NULL) + coord_cartesian(ylim = c(0,3e-5))) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
husler_reiss_low_b1_axis <- c(0, 1e-7)
husler_reiss_low_b2_axis <- c(0, 2.5e-6)
husler_reiss_low_b2_axis <- c(0, 3e-5)
ggsave("bma_update_deck/husler_reiss_low_fitted_k_not1.pdf",
       plot = husler_reiss_low_fitted_k_not1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)

husler_reiss_low_fitted_k_is1_all <- (low_b1_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_low_b1_axis)) + 
  (low_b2_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_low_b2_axis)) + 
  (low_b3_k_is1_fitted + ggtitle(NULL) + coord_cartesian(ylim = husler_reiss_low_b3_axis)) +
  plot_layout(guides = 'collect') & 
  theme(legend.background = element_rect(fill='transparent'),         panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("bma_update_deck/husler_reiss_low_fitted_k_is1_all.pdf",
       plot = husler_reiss_low_fitted_k_is1_all,
       dpi = 320,
       bg = "transparent",
       width = 16, height =7)
