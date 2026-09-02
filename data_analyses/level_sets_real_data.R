# =============================================================================
# Plots posterior level set (gauge function boundary) estimates on the real
# fire weather data for both stations. Overlays the data cloud, threshold,
# and posterior level set contours for all 6 gauge functions.
#
# Inputs:    fits_and_weights/post_params_joint/...qs (via extract_post_params_real_data.R)
#            data/raw/{data_type}_expo.qs
# Outputs:   figures/both_level_sets.png
# =============================================================================

library(tidyverse)
library(patchwork)
library(qs)
library(gaugeDependence)
library(ggnewscale)
library(grafify)
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

## functions to plot posterior unit level sets of all six gauge fits from the truncated and censored likelihoods ------
# Evaluate one fitted gauge g(w, 1-w) over a grid of angles, using the
# posterior-median dependence parameters (element 1 is alpha, so dropped).
post_level_set <- function(gauge_type, data_type, likelihood, w_sim) {
  params <- as.numeric(extract_post_params_radial(gauge_type, likelihood, data_type)[-1])
  gauge_fn <- get_gauge_function(gauge_type)
  return(gauge_fn(w_sim, 1 - w_sim, params) |> as_tibble() |> rename(values = V1) |> mutate(gauge = gauge_type, w_sim = w_sim))
}

# Evaluate all six gauges on a common angle grid and join to the observed
# data cloud (radii scaled by log(n) to match the level-set scale).
all_post_level_sets <- function(data_type, likelihood) {
  data <- qread(sprintf("data/raw/%s_expo.qs", data_type))$cloud_tib |> 
    select(r, w1, r0_w, high) |> 
    mutate(w_sim = seq(0, 1, length.out = n()),
           r = r/log(n()))
  gauge_library <- c("gauss", "inv_log", "rectangular", "logistic", "asym_log", "dirichlet")
  gw <- lapply(gauge_library, function(x) {
    post_level_set(gauge_type = x, data_type = data_type, likelihood = likelihood, w_sim = data$w_sim)
  }) |> 
    bind_rows() |>
    mutate(gauge = case_when(gauge == 'logistic' ~ 'Logistic',
                             gauge == 'gauss' ~ 'Gaussian',
                             gauge == 'inv_log' ~ 'Inv. logistic',
                             gauge == 'asym_log' ~ 'Asym. logistic',
                             gauge == 'dirichlet' ~ 'Dirichlet',
                             gauge == 'rectangular' ~ 'Rectangular',
                             .default = gauge),
           gauge = factor(gauge, levels = c("Gaussian", "Inv. logistic", "Rectangular", "Logistic", "Asym. logistic", "Dirichlet")))
  return(data |> 
           left_join(gw))
}

# Data cloud with the six fitted unit level sets overlaid as 1/g(w). Under the
# truncated likelihood the below-threshold points are hidden (transparent),
# since that fit only uses exceedances.
level_sets_plot <- function(data_type, likelihood) {
  level_sets_tib <- all_post_level_sets(data_type, likelihood)
  
  if(likelihood == "trunc") {
    pt_color <- c("transparent", "blue")
  } else {
    pt_color <- c("lightblue", "blue")
  }
  p <- level_sets_tib |> ggplot(aes(x = w1, y = r, color = high)) + 
    geom_point(size = rel(1.5)) +
    scale_color_manual(values = pt_color, guide = "none") +
    
    new_scale_color() +
    
    # fits of all 6 gauge functions
    geom_path(aes(x = w_sim, y = 1/values, group = gauge, color = gauge), linewidth = 1.2) +
    scale_color_grafify(palette = "r4") +
    
    theme_classic() +
    scale_x_continuous(expand = expansion(mult = c(0,0.03))) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.03))) + 
    theme(legend.position.inside = c(0.9, 0.9),
          panel.background = element_rect(fill='transparent', color = 'transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'),
          legend.background = element_rect(fill='transparent', color='transparent'),
          legend.text = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.4)),
          axis.title = element_text(size = rel(1.4))) +
    xlab(expression("W"["1"])) + ylab("R/log(n)") + labs(color = "")
  
  return(p)
}

both_lhoods_plot <- function(data_type) {
  p <- level_sets_plot(data_type, "trunc") + guide_area() + level_sets_plot(data_type, "cens") + 
    plot_layout(guides = 'collect') +
    plot_layout(widths = c(0.95, 0.5, 0.95)) &
    theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'))
  return(p)
}

## plot redstone and friendmtn together
friend_plots <- both_lhoods_plot("friendmtn")
redstone_plots <- both_lhoods_plot("redstone")

friend_plots / plot_spacer() / redstone_plots + plot_layout(heights = c(0.95, 0.05, 0.95))

ggsave("figures/both_level_sets.png",
       dpi = 320,
       bg = 'transparent',
       width = 12, height = 12)
knitr::plot_crop("figures/both_level_sets.png")
