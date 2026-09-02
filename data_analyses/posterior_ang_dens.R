# =============================================================================
# Plots the posterior angular density estimates (both star-shaped and mixture
# models) for each station and gauge, overlaid on the empirical angle histogram.
# Produces the angular density figures shown in the paper.
#
# Inputs:    fits_and_weights/post_params_joint/...qs (via extract_post_params_real_data.R)
#            data/raw/{data_type}_expo.qs
# Outputs:   figures/post_ang_dens.png, figures/post_ang_dens_fires.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(evd)
library(progressr)
library(RcppSimdJson)
library(gaugeDependence)
library(qs)
library(grafify)
library(patchwork)
source("extraction_scripts/extract_post_params_real_data.R")

gauge_functions <- list(
  gauss = gauss_gauge,
  inv_log = inv_log_gauge,
  rectangular = rectangular_gauge,
  logistic = logistic_gauge,
  asym_log = asym_log_gauge,
  dirichlet = dirichlet_gauge
)

# Grab gauge function by string
get_gauge_function <- function(type_str) {
  if (!type_str %in% names(gauge_functions)) {
    stop("Unknown gauge type: ", type_str)
  }
  return(gauge_functions[[type_str]])
}

# Monte Carlo estimate of the unit star-body volume, normalising the
# star-shaped angular density.
est_volume <- function(n = 100, pars, gauge_type) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  gauge_fn <- get_gauge_function(gauge_type)
  gx <- gauge_fn(grid[,1], grid[,2], pars)
  return(mean(gx <= 1))
}

# Star-shaped angular density f(w) = 1 / (g(w)^2 * 2 * vol(L)).
star_dens <- function(w1, pars, gauge_type) {
  w2 <- 1 - w1
  mc_star <- est_volume(n = 100, pars, gauge_type)
  gauge_fn <- get_gauge_function(gauge_type)
  gw <- gauge_fn(w1, w2, pars)
  return(1 / (gw^2 * 2 * mc_star))
}

# Angular mixture density evaluated at a single posterior draw.
mix_dens_one_iter <- function(w, params) {
  alphas <- as.numeric(params[grepl("alphastar", names(params))])
  betas <- as.numeric(params[grepl("betastar", names(params))])
  weights <- as.numeric(params[grepl("probs", names(params))])
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}

# Posterior mean angular mixture density, averaged over all draws.
mix_dens_all_iters <- function(data_type, w = NULL) {
  if(is.null(w)) w <- seq(0.001, 0.999, length.out = 1000)
  mix_params <- extract_post_params_ang_mix(data_type, FALSE)
  n_chains <- nrow(mix_params)
  dens_by_iter <- sapply(1:n_chains, function(i) mix_dens_one_iter(w, mix_params[i,]))
  return(rowMeans(dens_by_iter))
}

# Posterior star density for one gauge, using its posterior-median parameters.
star_by_gauge <- function(w, gauge, data_type) {
  post_star <- extract_post_params_ang_star(gauge, data_type, TRUE)
  return(star_dens(w, as.numeric(post_star), gauge) |> 
           as_tibble() |> rename(value = V1) |>
           mutate(method = gauge,
                  id = 1:n(),
                  w = w))
}

star_by_data <- function(data_type, w = NULL) {
  if(is.null(w)) w <- seq(0.001, 0.999, length.out = 1000)
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(lapply(gauge_library, function(x) star_by_gauge(w = w,
                                                         gauge = x,
                                                         data_type = data_type)) |> 
           bind_rows())
}

## diagnostic of angular dens fit ## ------
full_dens_by_data <- function(data_type) {
  w_sim <- seq(0.001, 0.999, length.out = 1000)
  dens <- star_by_data(data_type, w_sim)
  wts_star <- qread(sprintf("fits_and_weights/wts_joint_model/%s_cens_star.qs", data_type))
  wtd_dens <- suppressMessages(dens |> left_join(wts_star) |>
                                 mutate(stacking_preds = value * stacking,
                                        pseudo_boot = pseudobma_boot * value,
                                        pseudo_noboot = pseudobma_noboot * value) |>
                                 group_by(id, w) |>
                                 summarize(stacking_predictions = sum(stacking_preds),
                                           pseudobma_boot_preds = sum(pseudo_boot),
                                           pseudobma_noboot_preds = sum(pseudo_noboot),
                                           .groups = "drop")) |> 
    select(-id) |>
    pivot_longer(cols = -w, names_to = "method", values_to = "dens") |>
    mutate(method = case_when(grepl("stacking", method) ~ 'Stacking',
                              grepl("noboot", method) ~ 'Pseudo-BMA',
                              grepl("boot", method) ~ 'Pseudo-BMA+'),
           method = as.factor(method), 
           type = "Star-shaped",
           dataset = data_type)
  beta_mix_dens <- mix_dens_all_iters(data_type, w_sim) |> as_tibble() |>
    rename(dens = value) |>
    mutate(method = NA, type = "Mixture", w = w_sim, dataset = data_type)
  return(rbind(wtd_dens, beta_mix_dens))
}

datasets <- c("friendmtn", "redstone") |> as_tibble()
full_dens_angles <- pmap(datasets, function(value) {
  full_dens_by_data(value)
}) |> bind_rows()

obs_full_angles <- rbind(qread(sprintf("data/raw/%s_expo.qs", "friendmtn"))$W |> as_tibble() |> mutate(dataset = "friendmtn"),
                         qread(sprintf("data/raw/%s_expo.qs", "redstone"))$W |> as_tibble() |> mutate(dataset = "redstone"))


obs_full_angles |> 
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 0.025, aes(y = after_stat(density)), 
                 boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_line(data = ~subset(full_dens_angles, method %in% c("Pseudo-BMA+", NA)), aes(x = w, y = dens, color = type), linewidth = 1, alpha =0.85) +
  scale_color_grafify(palette = "all_grafify", ColSeq = FALSE) +
  facet_wrap(. ~ dataset, scales = "free") +
  labs(color = "Density type",
       x = expression("W"["1"]),
       y = "Density") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  theme_classic() +
  theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = 'transparent'),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.3)),
        legend.title = element_text(size = rel(1.3)),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.9),
        strip.text = element_blank(),
        panel.spacing.x = unit(0.95, "cm", data = NULL),
        panel.spacing.y = unit(0.95, "cm", data = NULL),
        legend.background = element_rect(fill = 'transparent', color = 'transparent'))

ggsave("figures/post_ang_dens_fires.png",
       dpi = 320,
       bg = 'transparent',
       width = 10,
       height = 5)
knitr::plot_crop("figures/post_ang_dens_fires.png")


## check distributional assumption of angles over quantile threshold ## ----
# Angles of all observations vs those above the fitted gauge threshold r0(w),
# used to check whether the angular distribution changes in the tail.
all_vs_top_angles_by_gauge <- function(gauge, data_type, lhood, q = 0.95) {
  data <- qread(sprintf("data/raw/%s_expo.qs", data_type))
  w <- data$W
  r <- data$R
  
  # get guage function
  gauge_fn <- get_gauge_function(gauge)
  
  # extract posterior params
  post_radial <- extract_post_params_radial(gauge = gauge, likelihood = lhood, data = data_type, summarize = TRUE)
  post_shape <- as.numeric(post_radial$alpha)
  post_dep <- as.numeric(post_radial[-1])
  # calculate posterior rate parameter
  gauge_vals <- gauge_fn(w, 1-w, post_dep)
  # return level quantiles
  r0w_post <- qgamma(0.95, shape = post_shape, rate = as.vector(gauge_vals))
  w_above <- tibble(angles = w[r > r0w_post], gauge = gauge)
  
  return(w_above)
}

all_angles_by_lhood_data <- function(data_type, lhood) {
  gauge_library <- c("gauss", "inv_log","rectangular", "logistic", "asym_log", "dirichlet")
  all_quant_angles <- lapply(gauge_library, function(x) {
    all_vs_top_angles_by_gauge(gauge = x,
                               data_type = data_type,
                               lhood = lhood)
  }) |> bind_rows() |>
    mutate(gauge = case_when(gauge == 'logistic' ~ 'Logistic',
                             gauge == 'gauss' ~ 'Gaussian',
                             gauge == 'inv_log' ~ 'Inv. logistic',
                             gauge == 'asym_log' ~ 'Asym. logistic',
                             gauge == 'dirichlet' ~ 'Dirichlet',
                             gauge == 'rectangular' ~ 'Rectangular',
                             .default = gauge),
           gauge = factor(gauge, levels = c("Gaussian", "Inv. logistic", "Rectangular", "Logistic", "Asym. logistic", "Dirichlet")))
  return(all_quant_angles)
}


quant_angles_plot <- function(data_type, lhood, ylims = c(NA, 2)) {
  data <- qread(sprintf("data/raw/%s_expo.qs", data_type))
  w <- data$W
  r <- data$R
  
  w_tib <- tibble(angles = w)
  
  quant_angles <- all_angles_by_lhood_data(data_type, lhood)
  
  w_tib |> 
    ggplot(aes(x = angles)) +
    geom_histogram(binwidth = 0.025, aes(y = after_stat(density)), 
                   boundary = 0, color = "gray45", fill = "lightgrey") +
    stat_density(data = quant_angles, 
                 aes(x = angles, color = gauge, group = gauge), 
                 geom = "line", position = "identity", linewidth = 1, alpha = 0.85) +
    scale_color_grafify(palette = "r4") +
    xlab(expression("W"["1"])) + ylab("Density") +
    scale_y_continuous(limits = ylims, expand = expansion(mult = c(0,0.01))) +
    scale_x_continuous(limits = c(NA, 1), expand = expansion(mult = c(0,0.04))) +
    labs(color = "Gauge function") +
    theme_classic() + 
    theme(panel.background = element_rect(fill='transparent', color='transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'),
          axis.text = element_text(size = rel(1.3)),
          axis.title = element_text(size = rel(1.5)),
          legend.text = element_text(size = rel(1.2)),
          legend.title = element_text(size = rel(1.2)),
          legend.position = "inside",
          legend.position.inside = c(0.9, 0.75),
          legend.background = element_rect(fill='transparent', color='transparent'))
}

quant_angles_plot("friendmtn", "cens", c(NA, 2.9)) + quant_angles_plot("friendmtn", "trunc", c(NA, 2.9))

both_lhoods_plot <- function(data_type, ylims = NULL) {
  p <- quant_angles_plot(data_type, "trunc", ylims) + plot_spacer() + (quant_angles_plot(data_type, "cens", ylims) + theme(legend.position = "none")) + 
    plot_layout(widths = c(1, 0.01, 1)) &
    theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'))
  return(p)
}

both_lhoods_plot("friendmtn", c(NA, 2.9)) / plot_spacer() / both_lhoods_plot("redstone", c(NA, 2.9)) + plot_layout(heights = c(0.95, 0.05, 0.95))

ggsave("figures/quant_assumpt_plot.pdf",
       bg = 'transparent',
       width = 10,
       height = 10,
       dpi = 320)
knitr::plot_crop("figures/quant_assumpt_plot.pdf")
