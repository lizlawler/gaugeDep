# =============================================================================
# Plots the distribution of angles W for observations above the threshold,
# comparing across dependence structures and levels. Used to illustrate
# how the angular distribution changes with dependence strength.
#
# Inputs:    data/{dep_type}/{dep_level}_{i}.json
# Outputs:   figures/angles_above_threhold.pdf
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(progressr) 
library(RcppSimdJson)
library(gaugeDependence)
library(qs)
library(grafify)
library(stringr)
library(scales)
library(patchwork)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

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

# For one dataset, return the angles W of (a) all observations and (b) those
# above the posterior gauge threshold under the censored and truncated
# radial fits. r0w = q-quantile of Gamma(alpha, g(W)) at each angle.
all_vs_top_angles <- function(gauge, level, datanum, q = 0.95) {
  data <- RcppSimdJson::fload(sprintf("data/%s/%s_%s.json", gauge, level, datanum))
  w <- data$W
  r <- data$R
  
  w_tib <- tibble(angles = w, type = "all")
  
  gauge_fn <- get_gauge_function(gauge)
  
  # with censored likelihood posterior params
  cens_params <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_cens_radial.qs", gauge, gauge, level))[datanum,]
  cens_dep <- cens_params[["dep"]]
  cens_alpha <- cens_params[["alpha"]]
  cens_gauge_vals <- gauge_fn(w, 1 - w, cens_dep)
  r0w_post_cens <- qgamma(q, shape = cens_alpha, rate = cens_gauge_vals)
  w_above_cens <- tibble(angles = w[r > r0w_post_cens], type = "Censored")
  
  # with truncated likelihood posterior params
  trunc_params <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_trunc_radial.qs", gauge, gauge, level))[datanum,]
  trunc_dep <- trunc_params[["dep"]]
  trunc_alpha <- trunc_params[["alpha"]]
  trunc_gauge_vals <- gauge_fn(w, 1 - w, trunc_dep)
  r0w_post_trunc <- qgamma(q, shape = trunc_alpha, rate = trunc_gauge_vals)
  w_above_trunc <- tibble(angles = w[r > r0w_post_trunc], type = "Truncated")

  angles_tib <- rbind(w_tib, w_above_cens, w_above_trunc) |>
    mutate(dep_type = sprintf("%s_%s", gauge, level))
  return(angles_tib)
}

# Representative dataset used for the illustration (one of the 200).
data_num <- 15
gauges <- c("gauss", "logistic")
levels <- c("low", "mid", "high")
all_combos <- expand_grid(gauges, levels)

all_angles <- pmap(all_combos, function(gauges, levels) {
  all_vs_top_angles(gauge = gauges, level = levels, datanum = data_num)
}) |> bind_rows() |>
  mutate(dep_type = factor(dep_type,
                           levels = c("gauss_low", "gauss_mid", "gauss_high",
                                      "logistic_low", "logistic_mid", "logistic_high")))

all_angles |> filter(type == "all") |> ggplot(aes(x = angles)) + 
  geom_histogram(data = ~subset(., dep_type == "gauss_low"), aes(y = after_stat(density)), binwidth = 0.04, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "gauss_mid"), aes(y = after_stat(density)), binwidth = 0.025, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "gauss_high"), aes(y = after_stat(density)), binwidth = 0.02, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "logistic_low"), aes(y = after_stat(density)), binwidth = 0.04, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "logistic_mid"), aes(y = after_stat(density)), binwidth = 0.025, boundary = 0, color = "gray45", fill = "lightgrey") +
  geom_histogram(data = ~subset(., dep_type == "logistic_high"), aes(y = after_stat(density)), binwidth = 0.01, boundary = 0, color = "gray45", fill = "lightgrey") +
  stat_density(data = ~subset(all_angles, type != "all"), aes(x = angles, group = type, color = type), geom = "line", position = "identity", linewidth = 1, alpha = 0.85) +
  scale_color_grafify(palette = "r4", ColSeq = FALSE) +
  facet_wrap(~ dep_type, scales = "free_y", axes = "all") +
  xlab(expression("W"["1"])) + ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  scale_x_continuous(limits = c(NA, 1), expand = expansion(mult = c(0,0.04))) +
  labs(color = "Likelihood") +
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

ggsave("figures/angles_above_threhold.pdf",
       bg = "transparent",
       dpi = 320,
       width = 12,
       height = 8)
knitr::plot_crop("figures/angles_above_threhold.pdf")
