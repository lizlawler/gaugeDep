# library(geometricMVE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(progressr)
library(patchwork)
source("gauge_functions_wrt_w.R")
source("run_wc_models.R")
options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

# reshape model fits and weights for easier future access (only need to do once) -------
# reshape stacking weights into usable format
# make_wts_df <- function(data_fit_type) {
#   weights_file <- list.files(path = "stacking_weights/",
#                              pattern = data_fit_type, 
#                              full.names = TRUE)
#   gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
#   temp <- readRDS(weights_file) |>
#     bind_rows() |> 
#     mutate(method = rep(gauge_library, 100)) |>
#     mutate(stacking = as.numeric(stacking),
#            pseudobma_boot = as.numeric(pseudobma_boot),
#            pseudobma_noboot = as.numeric(pseudobma_noboot)) |>
#     mutate(dataset = rep(1:100, times = rep(6, 100)))
#   return(temp)
# }
# 
# # reshape params
# reshape_params <- function(params_file, data_fit_type) {
#   model_name <- str_remove(str_remove(basename(params_file), "_all_iter_params.RDS"), paste0("_", data_fit_type))
#   params_list <- readRDS(params_file)
#   lapply(params_list, function(tib) {
#     temp <- tib |> select(-c(draw, dataset)) |> apply(2, median) |> as.list()
#     temp$dataset <- unique(tib$dataset)
#     return(temp)
#   }) |> 
#     bind_rows() |> 
#     mutate(median_params = pmap(pick(-last_col()), list)) |> 
#     select(last_col()-1, last_col()) |>
#     mutate(method = model_name)
# }
# 
# group_model_fits <- function(data_fit_type) {
#   all_fits <- list.files(path = "extracted_params/",
#                          pattern = data_fit_type, 
#                          full.names = TRUE)
#   temp <- lapply(all_fits, function(x) reshape_params(x, data_fit_type))
#   return(temp |> bind_rows())
# }
# 
# fits_and_weights <- function(data_fit_type) {
#   temp_wts <- make_wts_df(data_fit_type)
#   temp <- group_model_fits(data_fit_type) |> left_join(temp_wts)
#   dataset_lists <- split(temp, temp$dataset)
#   saveRDS(dataset_lists, 
#           file = paste0("fits_and_weights/", data_fit_type, ".RDS"))
#   print(paste0("Fits and weights have been reshaped and written to disk for ", data_fit_type))
# }
# 
# dep_types <- c("gauss", "logistic", "husler_reiss")
# dep_levels <- c("low", "mid", "high")
# lhood_type <- c("cens", "trunc")
# thresh_type <- c("marg", "ctau")
# 
# data_fit_combos <- expand_grid(dep_types, dep_levels, lhood_type, thresh_type) |> 
#   filter(!(dep_types == "husler_reiss" & thresh_type == "ctau"))
# 
# with_progress({
#   # Create a progress handler
#   p <- progressor(steps = nrow(data_fit_combos))
#   
#   # Apply the function using apply and update the progress bar
#   apply(data_fit_combos, 1, function(row) {
#     p()  # Update the progress bar
#     data_fit_type <- paste0(row["dep_types"], "_",
#                             row["dep_levels"], "_",
#                             row["lhood_type"], "_",
#                             row['thresh_type'])
#     fits_and_weights(data_fit_type)
#   })
# })


# simulate new bivariate data -------
sim_new_data_mixture <- function(dep_type, dep_level, data_num, 
                                 likelihood, threshold, 
                                 bma_method, nsim = 5000,
                                 fitted_gauge = F) {
  dataset <- RcppSimdJson::fload(paste0("data/", dep_type, "/", dep_level,"_", data_num, ".json"))
  w <- dataset$W
  r <- dataset$R
  r0w <- dataset$r0_w_ctau

  fit_wts <- readRDS(paste0("fits_and_weights/", dep_type, "_", 
                            dep_level, "_", 
                            likelihood, "_", 
                            threshold, ".RDS"))[[data_num]]
  mixture_prop <- sample(fit_wts$method, size = nsim, replace = TRUE, prob = fit_wts[[bma_method]]) |>
    as_tibble() |>
    group_by(value) |> 
    summarize(n = n())
  
  sim_data <- apply(mixture_prop, 1, function(row) {
    gauge_fcn <- get(paste0(row[["value"]], "_gauge"))
    pars <- fit_wts |> 
      filter(method == row[["value"]]) |> 
      select(median_params) |> 
      unlist(use.names = F)
    
    if(fitted_gauge) {
      # create threshold based on fitted gauge function
      gw_fitted <- gauge_fcn(w, pars[2:length(pars)])
      ctau_fitted <- quantile(gw_fitted * r, 0.95)
      r0w <- ctau_fitted/gw_fitted
    } else {
      # pull true gauge function threshold from data
      r0w <- data$r0_w_ctau
    }
    idx <- which(r > r0w)
    w <- w[idx]
    r0w <- r0w[idx]

    resample_idx <- sample(1:length(w), size = as.numeric(row[["n"]]), replace = T)
    r0w_resample <- r0w[resample_idx]
    w_resample <- w[resample_idx]
    rate_star <- gauge_fcn(w_resample, dep_par = pars[2:length(pars)])
    rstar <- qgamma(1 - runif(as.numeric(row[["n"]])) * pgamma(r0w_resample, shape = pars[1], 
                                             rate = rate_star, lower.tail = F),
                    shape = pars[1], rate = rate_star)
    xstar <- cbind(rstar * w_resample, rstar * (1 - w_resample))
    return(xstar |> as_tibble() |> rename(x1 = V1, x2 = V2))
  }) |> bind_rows()
  return(list(df = sim_data,
              dataset_num = data_num))
}

wc_sim_new_data <- function(dep_type, dep_level, data_num, nsim = 5000) {
  # read in data and recreate empirircal QR threshold
  dataset <- RcppSimdJson::fload(paste0("data/", dep_type, "/", dep_level,"_", data_num, ".json"))
  temp_qr <- geometricMVE::QR.2d(r = dataset$R, w = dataset$W, method = "empirical")
  r0w <- temp_qr$r0w   
  idx <- which(dataset$R > r0w)
  w <- dataset$W[idx]
  r0w <- r0w[idx]
  
  # pull chosen model for use in predictions
  temp <- split_wc_fits[[paste0(dep_type, ".", dep_level)]] |>
    filter(datasets == data_num)
  gauge_fcn <- get(paste0(temp$gauge_name, "_gauge"))
  pars <- unlist(temp$mle, use.names = F)
  
  resample_idx <- sample(1:length(w), size = nsim, replace = T)
  r0w_resample <- r0w[resample_idx]
  w_resample <- w[resample_idx]
  rate_star <- gauge_fcn(w_resample, dep_par = pars[2:length(pars)])
  rstar <- qgamma(1 - runif(nsim) * pgamma(r0w_resample, shape = pars[1], 
                                           rate = rate_star, lower.tail = F),
                  shape = pars[1], rate = rate_star)
  xstar <- cbind(rstar * w_resample, rstar * (1 - w_resample))
  return(xstar |> as_tibble() |> rename(x1 = V1, x2 = V2))
}

random_dataset <- sample(1:100, 1)
ls_new_data_cens_pseudoboot <- sim_new_data_mixture("husler_reiss", "mid", data_num = random_dataset, 
                                 likelihood = "cens", threshold = "marg",
                                 bma_method = "pseudobma_boot", nsim = 5000,fitted_gauge = T)
ls_new_data_trunc_pseudoboot <- sim_new_data_mixture("husler_reiss", "mid", data_num = random_dataset, 
                                         likelihood = "trunc", threshold = "marg",
                                         bma_method = "pseudobma_boot", nsim = 5000,fitted_gauge = T)
ls_new_data_cens_stacking <- sim_new_data_mixture("husler_reiss", "mid", data_num = random_dataset, 
                                                    likelihood = "cens", threshold = "marg",
                                                    bma_method = "stacking", nsim = 5000,fitted_gauge = T)
ls_new_data_trunc_stacking <- sim_new_data_mixture("husler_reiss", "mid", data_num = random_dataset, 
                                                     likelihood = "trunc", threshold = "marg",
                                                     bma_method = "stacking", nsim = 5000,fitted_gauge = T)
wc_new_data <- wc_sim_new_data("husler_reiss", "mid", random_dataset)

true_data <- RcppSimdJson::fload(paste0("data/husler_reiss/", "mid_", random_dataset, ".json"))
x1 <- true_data$R * true_data$W
x2 <- true_data$R * (1-true_data$W)
truth <- cbind(x1, x2) |> as_tibble()
cens_pseudoboot <- ggplot(truth, aes(x = x1, y = x2)) + geom_point(color = "darkgrey", alpha = 0.5, size = 0.3) + theme_classic() +
  geom_point(data = ls_new_data_cens_pseudoboot$df, aes(x = x1, y = x2), color = "blue", alpha = 0.5, size = 0.3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14)) + scale_y_continuous(expand = c(0,0), limits = c(0, 14)) +
  ggtitle("Censored, Pseudo-BMA+") + 
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA))

trunc_pseudoboot <- ggplot(truth, aes(x = x1, y = x2)) + geom_point(color = "darkgrey", alpha = 0.5, size = 0.3) + theme_classic() +
  geom_point(data = ls_new_data_trunc_pseudoboot$df, aes(x = x1, y = x2), color = "blue", alpha = 0.5, size = 0.3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14)) + scale_y_continuous(expand = c(0,0), limits = c(0, 14)) +
  ggtitle("Truncated, Pseudo-BMA+") + 
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA))

cens_stacking <- ggplot(truth, aes(x = x1, y = x2)) + geom_point(color = "darkgrey", alpha = 0.5, size = 0.3) + theme_classic() +
  geom_point(data = ls_new_data_cens_stacking$df, aes(x = x1, y = x2), color = "blue", alpha = 0.5, size = 0.3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14)) + scale_y_continuous(expand = c(0,0), limits = c(0, 14)) +
  ggtitle("Censored, Stacking") + 
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA))

trunc_stacking <- ggplot(truth, aes(x = x1, y = x2)) + geom_point(color = "darkgrey", alpha = 0.5, size = 0.3) + theme_classic() +
  geom_point(data = ls_new_data_trunc_stacking$df, aes(x = x1, y = x2), color = "blue", alpha = 0.5, size = 0.3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14)) + scale_y_continuous(expand = c(0,0), limits = c(0, 14)) + 
  ggtitle("Truncated, Stacking") + 
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA))

wc_plot <- ggplot(truth, aes(x = x1, y = x2)) + geom_point(color = "darkgrey", alpha = 0.5, size = 0.3) + theme_classic() +
  geom_point(data = wc_new_data, aes(x = x1, y = x2), color = "blue", alpha = 0.5, size = 0.3) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 14)) + scale_y_continuous(expand = c(0,0), limits = c(0, 14)) +
  ggtitle("W-C") + 
  theme(panel.background = element_rect(fill='transparent'), plot.background = element_rect(fill='transparent', color=NA))

# layout <- c(
#   area(t = 1, b = 2, l = 1, r = 2),
#   area(t = 1, b = 2, l = 3, r = 4),
#   area(t = 3, b = 4, l = 1, r = 2),
#   area(t = 3, b = 4, l = 3, r = 4),
#   area(t = 2, b = 3, l = 5, r = 6)
# )
# plot(layout)

all_sims_plot <- cens_pseudoboot + cens_stacking + trunc_pseudoboot + trunc_stacking + wc_plot + plot_layout(design = layout) 
all_sims_plot
ggsave("bma_update_deck/husler_reiss_mid_predictions.pdf",
       all_sims_plot,
       bg = 'transparent',
       width = 10,
       height = (20/3), dpi = 320)



