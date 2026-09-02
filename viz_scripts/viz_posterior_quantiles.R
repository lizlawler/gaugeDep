# =============================================================================
# Produces plots of posterior predictive quantiles of the gauge-function
# radial model against empirical quantiles, for calibration assessment
# across all 200 simulation datasets.
#
# Inputs:    fits_and_weights/post_params_joint/...qs
#            data/{dep_type}/{dep_level}_{i}.json
# Outputs:   figures/calibration_qs/{gauge}_{likelihood}_coverage.qs
#            figures/{gauge}_{dep_level}_coverage.pdf
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

viz_quant_reg <- function(gauge, level, datanum, q = 0.95) {
  data <- RcppSimdJson::fload(sprintf("data/%s/%s_%s.json", gauge, level, datanum))
  w <- data$W
  r <- data$R
  
  data_tib <- tibble(x = r * w, y = (1 - w) * r, dep_type = sprintf("%s_%s", gauge, level))
  
  gauge_fn <- get_gauge_function(gauge)
  
  # with censored likelihood posterior params
  cens_params <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_cens_radial.qs", gauge, gauge, level))[datanum,]
  cens_dep <- cens_params[["dep"]]
  cens_alpha <- cens_params[["alpha"]]
  
  cens_gauge_vals <- gauge_fn(w, 1 - w, cens_dep)
  r0w_post_cens <- qgamma(q, shape = cens_alpha, rate = cens_gauge_vals)
  
  # with truncated likelihood posterior params
  trunc_params <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_trunc_radial.qs", gauge, gauge, level))[datanum,]
  trunc_dep <- trunc_params[["dep"]]
  trunc_alpha <- trunc_params[["alpha"]]
  
  trunc_gauge_vals <- gauge_fn(w, 1 - w, trunc_dep)
  r0w_post_trunc <- qgamma(q, shape = trunc_alpha, rate = trunc_gauge_vals)
  
  thresh_tib <- tibble(x0 = r0w_post_cens * w, y0 = r0w_post_cens * (1 - w), thresh = "Censored") |> 
    rbind(tibble(x0 = r0w_post_trunc * w, y0 = r0w_post_trunc * (1 - w), thresh = "Truncated")) |> mutate(dep_type = sprintf("%s_%s", gauge, level))
  return(list(data_tib = data_tib, thresh_tib = thresh_tib))
}

data_num <- 15
gauss_high <- viz_quant_reg("gauss", "high", data_num)
gauss_mid <- viz_quant_reg("gauss", "mid", data_num)
gauss_low <- viz_quant_reg("gauss", "low", data_num)

logistic_high <- viz_quant_reg("logistic", "high", data_num)
logistic_mid <- viz_quant_reg("logistic", "mid", data_num)
logistic_low <- viz_quant_reg("logistic", "low", data_num)

all_data_types <- rbind(gauss_low$data_tib, gauss_mid$data_tib, gauss_high$data_tib, logistic_low$data_tib, logistic_mid$data_tib, logistic_high$data_tib) |>
  mutate(dep_type = factor(dep_type, levels = c("gauss_low", "gauss_mid", "gauss_high", "logistic_low", "logistic_mid", "logistic_high")))
all_thresh_types <- rbind(gauss_low$thresh_tib, gauss_mid$thresh_tib, gauss_high$thresh_tib, logistic_low$thresh_tib, logistic_mid$thresh_tib, logistic_high$thresh_tib) |>
  mutate(dep_type = factor(dep_type, levels = c("gauss_low", "gauss_mid", "gauss_high", "logistic_low", "logistic_mid", "logistic_high")))

ggplot(all_data_types, aes(x = x, y = y)) + geom_point(alpha = 0.5, size = 1, color = "black") +
  geom_point(data = all_thresh_types, aes(x = x0, y = y0, color = thresh, group = thresh), size = 0.5) +
  scale_color_grafify(palette = "r4", ColSeq = FALSE) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"])) +
  scale_y_continuous(limits = c(NA, 10), expand = expansion(mult = c(0,0.03))) +
  scale_x_continuous(limits = c(NA, 10), expand = expansion(mult = c(0,0.04))) +
  facet_wrap(. ~ dep_type, axes = "all") +
  theme_classic() + 
  guides(color = guide_legend(override.aes = list(size = 3) ) ) +
  labs(color = "Likelihood") +
  theme(panel.background = element_rect(fill='transparent', color='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.5)),
        panel.spacing.x = unit(0.95, "cm", data = NULL),
        panel.spacing.y = unit(0.95, "cm", data = NULL),
        strip.text.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.925, 0.075),
        legend.text = element_text(size = rel(1.3)),
        legend.title = element_text(size = rel(1.3)),
        legend.background = element_rect(fill='transparent', color='transparent'))

ggsave("figures/quantile_threshold.png",
       dpi = 320,
       bg = "transparent",
       width = 12,
       height = 8)
knitr::plot_crop("figures/quantile_threshold.png")
