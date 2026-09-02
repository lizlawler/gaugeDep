# =============================================================================
# Visualises the posterior angular density estimates across all 6 gauge
# functions for the simulation study datasets, comparing star-shaped and
# mixture model fits. Produces density overlay plots for representative datasets.
#
# Inputs:    fits_and_weights/post_params_joint/...qs
#            data/{dep_type}/{dep_level}_{i}.json
# Outputs:   figures/viz_ang_dens.png
# =============================================================================

library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(cmdstanr)
library(patchwork)
library(posterior)
library(grafify)
library(gaugeDependence)
library(qs)

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

est_volume <- function(n = 100, pars = 0.5, gauge_type) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  return(est_star_vol(grid_x = as.matrix(grid), pars = pars, ang_gauge_type = gauge_type))
}

# Star-shaped angular density f(w) = 1 / (g(w)^2 * 2 * vol(L)). The star-body
# volume vol(L) is estimated by Monte Carlo, except for the logistic gauge
# where it equals the dependence parameter directly.
star_dens <- function(w1, pars, gauge_type) {
  w2 <- 1 - w1
  mc_star <- if(gauge_type!= "logistic") {
    est_volume(n = 100, pars, gauge_type) 
  } else {
    pars
  }
  gauge_fn <- get_gauge_function(gauge_type)
  gw <- gauge_fn(w1, w2, pars)
  return(1 / (gw^2 * 2 * mc_star))
}

# Angular mixture density: weighted sum of Beta components from the NIMBLE fit.
mix_dens <- function(w, mean_params) {
  alphas <- as.numeric(mean_params[grepl("alphastar", names(mean_params))])
  betas <- as.numeric(mean_params[grepl("betastar", names(mean_params))])
  weights <- as.numeric(mean_params[grepl("probs", names(mean_params))])
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}

# data_num <- 15
# read in MCMC samples of different fits
calc_ang_density <- function(gauge, level, datanum) {
  data <- RcppSimdJson::fload(sprintf("data/%s/%s_%s.json", gauge, level, datanum))
  w <- data$W
  
  mix_params <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_ang_mix.qs", gauge, level))[datanum,]
  star_params <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_ang_star.qs", gauge, level, gauge))$dep[datanum]
  
  w_sim <- seq(0, 1, length.out = length(w))
  mix_curve <- as.numeric(mix_dens(w_sim, mix_params))
  star_curve <- as.numeric(star_dens(w_sim, star_params, gauge_type = gauge))
  
  dens_tib <- cbind(w, w_sim) |> as_tibble() |>
    mutate("Mixture" = mix_curve,
           "Star-shaped" = star_curve) |>
    pivot_longer(cols = -c(w, w_sim), names_to = "Density", values_to = "dens_val") |>
    mutate(Density = as.factor(Density),
           dep_type = sprintf("%s_%s", gauge, level))

  return(dens_tib)
}

# Representative dataset used for the illustration (one of the 200).
data_num <- 62

gauges <- c("gauss", "logistic")
levels <- c("low", "mid", "high")
all_combos <- expand_grid(gauges, levels)

all_dens <- pmap(all_combos, function(gauges, levels) {
  calc_ang_density(gauge = gauges, level = levels, datanum = data_num)
}) |> bind_rows() |>
  mutate(dep_type = factor(dep_type,
                           levels = c("gauss_low", "gauss_mid", "gauss_high",
                                      "logistic_low", "logistic_mid", "logistic_high")))


full_plot <- ggplot(all_dens, aes(x = w)) + 
  geom_histogram(data = ~subset(., dep_type == "gauss_low"), aes(y = after_stat(density)), binwidth = 0.04, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "gauss_mid"), aes(y = after_stat(density)), binwidth = 0.025, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "gauss_high"), aes(y = after_stat(density)), binwidth = 0.02, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "logistic_low"), aes(y = after_stat(density)), binwidth = 0.04, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "logistic_mid"), aes(y = after_stat(density)), binwidth = 0.025, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "logistic_high"), aes(y = after_stat(density)), binwidth = 0.01, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_line(aes(x = w_sim, y = dens_val, color = Density), linewidth = 1, alpha = 0.9) +
  scale_color_grafify(palette = "all_grafify", ColSeq = FALSE) +
  facet_wrap(~ dep_type, scales = "free_y", axes = "all") +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  scale_x_continuous(limits = c(NA, 1), expand = expansion(mult = c(0,0.04))) +
  labs(color = "Density type") +
  theme_classic() + 
  theme(panel.background = element_rect(fill='transparent', color='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.3)),
        legend.title = element_text(size = rel(1.3)),
        strip.text = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.935, 0.35),
        panel.spacing.x = unit(0.95, "cm", data = NULL),
        panel.spacing.y = unit(0.95, "cm", data = NULL),
        legend.background = element_rect(fill='transparent', color='transparent'))


ggsave("figures/viz_ang_dens.png",
       plot = full_plot,
       dpi = 320,
       bg = "transparent",
       width = 12, height = 8)
knitr::plot_crop("figures/viz_ang_dens.png")



