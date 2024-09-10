library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(cmdstanr)
library(patchwork)
library(posterior)
# library(dirichletprocess)
source("gauge_functions_wrt_x.R")

inv_logit <- function(x) {
  return(1 / (1 + exp(-x)))
}

approx_indicator <- function(x, k) {
  temp_arg <- k * (1.0 - x)
  return(inv_logit(k * (1.0 - x)))
}

est_volume <- function(n = 100, pars = 0.5, gauge, approx = FALSE) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  gx <- gauge_fcn(grid[,1], grid[,2], dep_par = pars)
  if(approx) {
    return(mean(approx_indicator(gx, 15)))
  } else {
    return(mean(gx <= 1))
  }
}

# density functions ---------
# L1, pseudo angles
dens_l1_norm <- function(w1, w2, gauge, par_val, approx = FALSE) {
  if(gauge == "logistic") {
    mc_vol <- par_val
  } else if(approx) {
    mc_vol <- est_volume(n = 100, par_val, gauge, approx = TRUE)
  } else {
    mc_vol <- est_volume(n = 100, par_val, gauge, approx = FALSE)
  }
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  gw <- gauge_fcn(w1, w2, par_val)
  return(1 / (gw^2 * 2 * mc_vol))
}

bernstein_dens <- function(w, weights) {
  k <- length(weights)
  pdf <- 0.0
  for(j in 1:k) {
    pdf = pdf + weights[j] * dbeta(w, shape1 = j, shape2 = (k - j + 1))
  }
  return(pdf)
}

dpmm_density <- function(w, weights, params) {
  k <- length(weights)
  pdf <- 0.0 
  for(j in 1:k) {
    alpha_post <- params_list[[1]][j] * params_list[[2]][j]
    beta_post <- (1-params_list[[1]][j]) * params_list[[2]][j]
    pdf = pdf + weights[j] * dbeta(w, shape1 = alpha_post, shape2 = beta_post)
  }
  return(pdf)
}
# logistic_high_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "high_1_bern_\\d{1}", full.names = TRUE)
# post_param <- apply(read_cmdstan_csv(logistic_high_1_files, variables = "weights")$post_warmup_draws, 3, median)
# 
# logistic_low_11_data <- RcppSimdJson::fload("data/angular/logistic/low_11.json")
# w1 <- logistic_low_1_data$w1
# hist(logistic_high_11_data$w1, freq = FALSE, breaks = 70)
# test_dens <- bernstein_dens(w1, post_param)
# points(w1, test_dens)

## gaussian angular density visualizations ---------
# low dependence, volume of star shaped set
gauss_low_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "low_1_\\d{1}", full.names = TRUE)
post_param <- median(read_cmdstan_csv(gauss_low_1_files, variables = "dep")$post_warmup_draws |> c())

gauss_low_1_data <- RcppSimdJson::fload("data/angular/gauss/low_1.json")
w1 <- gauss_low_1_data$w1

gauss_low_1_bern_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "low_1_bern", full.names = TRUE)
post_param_bern <- apply(read_cmdstan_csv(gauss_low_1_bern_files, variables = "weights")$post_warmup_draws, 3, median)

gauss_low_1_dpmm <- readRDS("dpmm_fits/dp_gauss_low.RDS")
params_list <- gauss_low_1_dpmm$clusterParameters
dpmm_weights <- gauss_low_1_dpmm$weights

# hist(w1, freq = FALSE, ylim = c(0,3))
# points(w1, dens_l1_norm(w1, (1-w1), "gauss", post_param))
low_gauss_tib <- tibble(angle = w1, 
                        est_dens = dens_l1_norm(w1, (1-w1), "gauss", post_param, approx = TRUE),
                        og_dens = dens_l1_norm(w1, (1-w1), "gauss", 0.1),
                        bern_dens = bernstein_dens(w1, post_param_bern),
                        dpmm_dens = dpmm_density(w1, dpmm_weights, params_list))
plot_low_gauss <- low_gauss_tib |> ggplot(aes(x = angle)) + 
  geom_histogram(binwidth = 0.05, aes(y = after_stat(density)), boundary = 0, color = "black", fill = "darkgrey") +
  geom_line(aes(x = angle, y = og_dens), color = "darkgoldenrod1", linewidth = 1) +
  geom_line(aes(x = angle, y = est_dens), color = "blueviolet", linewidth = 1) +
  geom_line(aes(x = angle, y = bern_dens), color = "salmon", linewidth = 1) +
  geom_line(aes(x = angle, y = dpmm_dens), color = "cornflowerblue", linewidth = 1) +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  theme_classic()
plot_low_gauss

# mid dependence
gauss_mid_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "mid_1_\\d{1}", full.names = TRUE)
post_param <- median(read_cmdstan_csv(gauss_mid_1_files, variables = "dep")$post_warmup_draws |> c())

gauss_mid_1_data <- RcppSimdJson::fload("data/angular/gauss/mid_1.json")
w1 <- gauss_mid_1_data$w1

gauss_mid_1_bern_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "mid_1_bern", full.names = TRUE)
post_param_bern <- apply(read_cmdstan_csv(gauss_mid_1_bern_files, variables = "weights")$post_warmup_draws, 3, median)

gauss_mid_1_dpmm <- readRDS("dpmm_fits/dp_gauss_mid.RDS")
params_list <- gauss_mid_1_dpmm$clusterParameters
dpmm_weights <- gauss_mid_1_dpmm$weights

# hist(w1, freq = FALSE, ylim = c(0,3))
# points(w1, dens_l1_norm(w1, (1-w1), "gauss", post_param))
mid_gauss_tib <- tibble(angle = w1, 
                        est_dens = dens_l1_norm(w1, (1-w1), "gauss", post_param, approx = TRUE),
                        og_dens = dens_l1_norm(w1, (1-w1), "gauss", 0.5),
                        bern_dens = bernstein_dens(w1, post_param_bern),
                        dpmm_dens = dpmm_density(w1, dpmm_weights, params_list))
plot_mid_gauss <- mid_gauss_tib |> ggplot(aes(x = angle)) + 
  geom_histogram(binwidth = 0.05, aes(y = after_stat(density)), boundary = 0, color = "black", fill = "darkgrey") +
  geom_line(aes(x = angle, y = og_dens), color = "darkgoldenrod1", linewidth = 1) +
  geom_line(aes(x = angle, y = est_dens), color = "blueviolet", linewidth = 1) +
  geom_line(aes(x = angle, y = bern_dens), color = "salmon", linewidth = 1) +
  geom_line(aes(x = angle, y = dpmm_dens), color = "cornflowerblue", linewidth = 1) +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  theme_classic()
plot_mid_gauss

# high dependence
gauss_high_1_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "high_1_\\d{1}", full.names = TRUE)
post_param <- median(read_cmdstan_csv(gauss_high_1_files, variables = "dep")$post_warmup_draws |> c())

gauss_high_1_data <- RcppSimdJson::fload("data/angular/gauss/high_1.json")
w1 <- gauss_high_1_data$w1

gauss_high_1_bern_files <- list.files(path = "stan/angular/csv_fits/gauss", pattern = "high_1_bern_\\d{1}", full.names = TRUE)
post_param_bern <- apply(read_cmdstan_csv(gauss_high_1_bern_files, variables = "weights")$post_warmup_draws, 3, median)

gauss_high_1_dpmm <- readRDS("dpmm_fits/dp_gauss_high.RDS")
params_list <- gauss_high_1_dpmm$clusterParameters
dpmm_weights <- gauss_high_1_dpmm$weights

# hist(w1, freq = FALSE, ylim = c(0,3))
# points(w1, dens_l1_norm(w1, (1-w1), "gauss", post_param))
high_gauss_tib <- tibble(angle = w1, 
                        est_dens = dens_l1_norm(w1, (1-w1), "gauss", post_param, approx = TRUE),
                        og_dens = dens_l1_norm(w1, (1-w1), "gauss", 0.9),
                        bern_dens = bernstein_dens(w1, post_param_bern),
                        dpmm_dens = dpmm_density(w1, dpmm_weights, params_list))
plot_high_gauss <- high_gauss_tib |> ggplot(aes(x = angle)) + 
  geom_histogram(binwidth = 0.05, aes(y = after_stat(density)), boundary = 0, color = "black", fill = "darkgrey") +
  geom_line(aes(x = angle, y = og_dens), color = "darkgoldenrod1", linewidth = 1) +
  geom_line(aes(x = angle, y = est_dens), color = "blueviolet", linewidth = 1) +
  geom_line(aes(x = angle, y = bern_dens), color = "salmon", linewidth = 1) +
  geom_line(aes(x = angle, y = dpmm_dens), color = "cornflowerblue", linewidth = 1) +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  theme_classic()
plot_high_gauss

row_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = "Gaussian",
                                                       face = "bold",
                                                       size = 10))
row_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = "Logistic",
                                                       face = "bold",
                                                       size = 10))
blank_space <- wrap_elements(panel = ggpubr::text_grob(label = '',
                                                       size = 5))

all_gauss_plots <- (row_label_1 | plot_low_gauss | plot_mid_gauss | plot_high_gauss) + plot_layout(widths = c(0.15, .95, .95, .95))

## logistic - angular density ---------
# low dependence
logistic_low_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "low_1_\\d{1}", full.names = TRUE)
post_param <- median(read_cmdstan_csv(logistic_low_1_files, variables = "dep")$post_warmup_draws |> c())

logistic_low_1_data <- RcppSimdJson::fload("data/angular/logistic/low_1.json")
w1 <- logistic_low_1_data$w1

logistic_low_1_bern_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "low_1_bern_\\d{1}", full.names = TRUE)
post_param_bern <- apply(read_cmdstan_csv(logistic_low_1_bern_files, variables = "weights")$post_warmup_draws, 3, median)

logistic_low_1_dpmm <- readRDS("dpmm_fits/dp_logistic_low.RDS")
params_list <- logistic_low_1_dpmm$clusterParameters
dpmm_weights <- logistic_low_1_dpmm$weights

low_logistic_tib <- tibble(angle = w1, 
                        est_dens = dens_l1_norm(w1, (1-w1), "logistic", post_param),
                        og_dens = dens_l1_norm(w1, (1-w1), "logistic", 0.9),
                        bern_dens = bernstein_dens(w1, post_param_bern),
                        dpmm_dens = dpmm_density(w1, dpmm_weights, params_list))
plot_low_logistic <- low_logistic_tib |> ggplot(aes(x = angle)) + 
  geom_histogram(binwidth = 0.05, aes(y = after_stat(density)), boundary = 0, color = "black", fill = "darkgrey") +
  geom_line(aes(x = angle, y = og_dens), color = "darkgoldenrod1", linewidth = 1) +
  geom_line(aes(x = angle, y = est_dens), color = "blueviolet", linewidth = 1) +
  geom_line(aes(x = angle, y = bern_dens), color = "salmon", linewidth = 1) +
  geom_line(aes(x = angle, y = dpmm_dens), color = "cornflowerblue", linewidth = 1) +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  theme_classic()
plot_low_logistic 

# mid dependence
logistic_mid_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "mid_1_\\d{1}", full.names = TRUE)
post_param <- median(read_cmdstan_csv(logistic_mid_1_files, variables = "dep")$post_warmup_draws |> c())

logistic_mid_1_data <- RcppSimdJson::fload("data/angular/logistic/mid_1.json")
w1 <- logistic_mid_1_data$w1

logistic_mid_1_bern_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "mid_1_bern_\\d{1}", full.names = TRUE)
post_param_bern <- apply(read_cmdstan_csv(logistic_mid_1_bern_files, variables = "weights")$post_warmup_draws, 3, median)

logistic_mid_1_dpmm <- readRDS("dpmm_fits/dp_logistic_mid.RDS")
params_list <- logistic_mid_1_dpmm$clusterParameters
dpmm_weights <- logistic_mid_1_dpmm$weights

mid_logistic_tib <- tibble(angle = w1, 
                           est_dens = dens_l1_norm(w1, (1-w1), "logistic", post_param),
                           og_dens = dens_l1_norm(w1, (1-w1), "logistic", 0.5),
                           bern_dens = bernstein_dens(w1, post_param_bern),
                           dpmm_dens = dpmm_density(w1, dpmm_weights, params_list))
plot_mid_logistic <- mid_logistic_tib |> ggplot(aes(x = angle)) + 
  geom_histogram(binwidth = 0.025, aes(y = after_stat(density)), boundary = 0, color = "black", fill = "darkgrey") +
  geom_line(aes(x = angle, y = og_dens), color = "darkgoldenrod1", linewidth = 1) +
  geom_line(aes(x = angle, y = est_dens), color = "blueviolet", linewidth = 1) +
  geom_line(aes(x = angle, y = bern_dens), color = "salmon", linewidth = 1) +
  geom_line(aes(x = angle, y = dpmm_dens), color = "cornflowerblue", linewidth = 1) +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  theme_classic()
plot_mid_logistic


# high dependence
logistic_high_1_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "high_1_mix_betas_take3_\\d{1}", full.names = TRUE)
# logistic_high_1_files <- logistic_high_1_files[!grepl("gauss", logistic_high_1_files)]
post_param <- read_cmdstan_csv(logistic_high_1_files,variables = c("alpha", "beta", "weights"))$post_warmup_draws
param_vec <-  post_param |> as_draws_df() |> select(-c(".chain", ".iteration", ".draw")) |> apply(2, median)

mix_dens <- function(w, alphas, betas, weights) {
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}
alphas <- param_vec[grepl("alpha", names(param_vec))]
betas <- param_vec[grepl("beta", names(param_vec))]
weights <- param_vec[grepl("weight", names(param_vec))]
logistic_high_1_data <- RcppSimdJson::fload("data/angular/logistic/high_1.json")
w1 <- logistic_high_1_data$w1
hist(w1, freq = FALSE, breaks = 60)
points(w1, mix_dens(w1, alphas, betas, weights), col = "blue")
logistic_high_1_bern_files <- list.files(path = "stan/angular/csv_fits/logistic", pattern = "high_1_bern_v2", full.names = TRUE)
post_param_bern <- apply(read_cmdstan_csv(logistic_high_1_bern_files, variables = "weights")$post_warmup_draws, 3, median)

logistic_high_1_dpmm <- readRDS("dpmm_fits/dp_logistic_high.RDS")
params_list <- logistic_high_1_dpmm$clusterParameters
dpmm_weights <- logistic_high_1_dpmm$weights

high_logistic_tib <- tibble(angle = w1, 
                           est_dens = dens_l1_norm(w1, (1-w1), "logistic", post_param),
                           og_dens = dens_l1_norm(w1, (1-w1), "logistic", 0.1),
                           bern_dens = bernstein_dens(w1, post_param_bern),
                           dpmm_dens = dpmm_density(w1, dpmm_weights, params_list))
plot_high_logistic <- high_logistic_tib |> ggplot(aes(x = angle)) + 
  geom_histogram(binwidth = 0.01, aes(y = after_stat(density)), boundary = 0, color = "black", fill = "darkgrey") +
  geom_line(aes(x = angle, y = og_dens), color = "darkgoldenrod1", linewidth = 1) +
  geom_line(aes(x = angle, y = est_dens), color = "blueviolet", linewidth = 1) +
  geom_line(aes(x = angle, y = bern_dens), color = "salmon", linewidth = 1) +
  geom_line(aes(x = angle, y = dpmm_dens), color = "cornflowerblue", linewidth = 1) +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  theme_classic()
plot_high_logistic

all_logistic_plots <- (row_label_2 | plot_low_logistic | plot_mid_logistic | plot_high_logistic) + plot_layout(widths = c(0.15, .95, .95, .95))

all_plots <- (all_gauss_plots / all_logistic_plots) + 
  plot_layout(heights = c(0.85, 0.85)) &
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("figures/angular_densities_bern_dpmm.pdf",
       plot = all_plots,
       dpi = 320,
       bg = "transparent",
       width = 15, height = 8)
