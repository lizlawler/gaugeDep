# =============================================================================
# Plots the unit-level-set boundaries {x : g(x) = 1} for all 6 gauge
# functions at a representative dependence parameter, illustrating how
# each gauge defines a different shape of extremal dependence.
#
# Inputs:    (uses gauge functions from the gaugeDependence package)
# Outputs:   figures/level_sets_gauges_labeled.png
# =============================================================================

library(tidyverse)
library(gaugeDependence)

w_gw_tib <- tibble(w = seq(0, 1, length.out = 500)) |> 
  mutate("Gaussian" = as.numeric(gauss_gauge(w, 1-w, 0.5)),
         "Logistic" = as.numeric(logistic_gauge(w, 1-w, 0.5)),
         "Inverted logistic" = as.numeric(inv_log_gauge(w, 1-w, 0.5)),
         "Rectangular" = as.numeric(rectangular_gauge(w, 1-w, 0.5)),
         "Asymmetric logistic" = as.numeric(asym_log_gauge(w, 1-w, 0.5)),
         "Dirichlet" = as.numeric(dirichlet_gauge(w, 1-w, c(0.25,3)))) |>
  pivot_longer(cols = -w, names_to = "gw", values_to = "value") |>
  mutate(gw = factor(gw, levels = c("Gaussian", "Inverted logistic", "Rectangular", "Logistic", "Asymmetric logistic", "Dirichlet")))

w_gw_tib |> ggplot(aes(x = w/value, y = (1-w)/value)) + 
  geom_path(color = "blue", linewidth = 1) + 
  geom_vline(xintercept = 1, color = "grey34", linetype = "dashed", linewidth = 0.75) +
  geom_hline(yintercept = 1, color = "grey34", linetype = "dashed", linewidth = 0.75) +
  facet_wrap(. ~ gw, axes= "all", axis.labels = "margins") +
  theme_classic() + coord_fixed() +
  scale_x_continuous(expand = expansion(mult = c(0,0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.02))) +
  # xlab("X") + ylab("Y") +
  theme(panel.background = element_rect(fill='transparent', color='black'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_blank(),
        # axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        # strip.text = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        strip.background = element_rect(fill='transparent', color='black'),
        panel.spacing.x = unit(0.95, "cm", data = NULL),
        panel.spacing.y = unit(0.95, "cm", data = NULL))

ggsave("figures/level_sets_gauges_labeled.png",
       bg = "transparent",
       dpi = 320,
       width = 12,
       height = 8)
knitr::plot_crop("figures/level_sets_gauges_labeled.png")
