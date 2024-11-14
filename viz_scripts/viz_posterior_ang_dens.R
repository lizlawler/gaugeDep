library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(cmdstanr)
library(patchwork)
library(posterior)
# library(dirichletprocess)
source("gauge_functions_wrt_x.R")
source("gauge_functions_wrt_w.R")

est_volume <- function(n = 100, pars = 0.5) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  # gauge_fcn <- get(paste0(gauge, "_gauge"))
  gx <- gauss_gauge_wrt_x(grid[,1], grid[,2], dep_par = pars)
  return(mean(gx <= 1))
}

# density functions ---------
# L1, pseudo angles
dens_l1_norm <- function(w1, par_val) {
  # if(gauge == "logistic") {
  #   mc_vol <- par_val
  # } else {
  mc_vol <- est_volume(n = 100, par_val)
  # }
  # gauge_fcn <- get(paste0(gauge, "_gauge"))
  gw <- gauss_gauge(w1, par_val)
  return(1 / (gw^2 * 2 * mc_vol))
}

mix_dens <- function(w, chain_of_params) {
  alphas <- as.numeric(chain_of_params[grepl("alpha", names(chain_of_params))])
  betas <- as.numeric(chain_of_params[grepl("beta", names(chain_of_params))])
  weights <- as.numeric(chain_of_params[grepl("weight", names(chain_of_params))])
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}

visualize_fit_density <- function(gauge, level, datanum, binwidth = 0.05) {
  data <- RcppSimdJson::fload(paste0("data/", gauge, "/", level, "_", datanum, ".json"))
  W <- data$W
  beta_mixture_files <- list.files(path = paste0("stan/radial_angular/csv_fits/", gauge, "/"), 
                          pattern = paste0(level, "_", datanum, "_\\d{1}.csv"), full.names = TRUE)
  beta_mixture_params <-  as_cmdstan_fit(beta_mixture_files) |> as_draws_df() |> 
    select(any_of(contains(c("weights", "alpha","beta", "dep", ".draw", ".chain")))) |>
    group_by(.chain) |> 
    summarize(across(everything(), median)) |> 
    select(-c(.draw, alpha, dep))
  beta_mixture_params_list <- split(beta_mixture_params, beta_mixture_params$.chain)
  starshaped_fit <- readRDS(paste0("mcmc_samples/", gauge, "/", level, "_", datanum, ".rds"))
  starshaped_params <- lapply(lapply(starshaped_fit, function(x) x$dep_w[10000:25000]), median)
  dens_tib <- as_tibble(W) |> rename(w = value) |>
    mutate(beta_chain1 = mix_dens(w, beta_mixture_params_list[[1]]),
           beta_chain2 = mix_dens(w, beta_mixture_params_list[[2]]),
           beta_chain3 = mix_dens(w, beta_mixture_params_list[[3]]),
           star_chain1 = dens_l1_norm(w, starshaped_params[[1]]),
           star_chain2 = dens_l1_norm(w, starshaped_params[[2]]),
           star_chain3 = dens_l1_norm(w, starshaped_params[[3]])) |>
    pivot_longer(cols = -w, names_to = "dens_type", values_to = "density") |>
    mutate(dens_type = as.factor(dens_type))
  dens_tib |> ggplot(aes(x = w)) + 
    geom_histogram(binwidth = binwidth, aes(y = after_stat(density)), 
                   boundary = 0, color = "black", fill = "lightgrey") +
    geom_line(aes(x = w, y = density, color = dens_type), linewidth = 1, alpha = 0.8) +
    xlab(expression("W"["1"])) + ylab("Density") +
    scale_y_continuous(expand = c(0,0.01)) +
    scale_x_continuous(expand = c(0.01,0.01)) +
    theme_classic()
}

visualize_fit_density("gauss", "mid", 15)
visualize_fit_density("gauss", "low", 15)
visualize_fit_density("gauss", "high", 15)

visualize_fit_density("logistic", "mid", 15)
visualize_fit_density("gauss", "low", 15)
visualize_fit_density("gauss", "high", 15)

# row_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = "Gaussian",
#                                                        face = "bold",
#                                                        size = 10))
# row_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = "Logistic",
#                                                        face = "bold",
#                                                        size = 10))
# blank_space <- wrap_elements(panel = ggpubr::text_grob(label = '',
#                                                        size = 5))
# 
# all_gauss_plots <- (row_label_1 | plot_low_gauss | plot_mid_gauss | plot_high_gauss) + plot_layout(widths = c(0.15, .95, .95, .95))
# 
# all_logistic_plots <- (row_label_2 | plot_low_logistic | plot_mid_logistic | plot_high_logistic) + plot_layout(widths = c(0.15, .95, .95, .95))
# 
# all_plots <- (all_gauss_plots / all_logistic_plots) + 
#   plot_layout(heights = c(0.85, 0.85)) &
#   theme(panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent'))
# 
# ggsave("figures/angular_densities_bern_dpmm.pdf",
#        plot = all_plots,
#        dpi = 320,
#        bg = "transparent",
#        width = 15, height = 8)
