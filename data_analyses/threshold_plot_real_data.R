# =============================================================================
# Produces the threshold selection and data cloud visualisation for the real
# fire weather data, showing both Euclidean (scatter plot) and polar (R vs W)
# coordinates with above-threshold points highlighted.
#
# Inputs:    data/raw/redstone_expo.qs
# Outputs:   figures/redstone_threshold_plots.pdf
# =============================================================================

library(evd)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(qs)

# Point colours: below threshold, above threshold, and the threshold curve.
cols <- c("lightblue", "blue", "red")
data <- qread("data/raw/redstone_expo.qs")
cloud_tib <- data$cloud_tib
n <- nrow(cloud_tib)

# Euclidean view: ERC vs FWI, scaled by log(n) so the threshold curve is
# on a comparable scale across sample sizes.
(data_euc_plot <- cloud_tib |>
  ggplot(aes(x / log(n), y / log(n), color = high)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(x = x_lb / log(n), y = y_lb / log(n), color = "red"), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_x_continuous(limits = c(0, 1.25), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(0, 1.25), expand = expansion(mult = c(0, 0.01))) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.position = "none"
  ) +
  xlab(expression("ERC" / "log(n)")) +
  ylab(expression("FWI" / "log(n)")))

# Pseudo-polar view: radius R against angle W1, with the same threshold
# curve expressed as r0(w).
(data_polar_plot <- cloud_tib |>
  ggplot(aes(w1, r / log(n), color = high)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(x = w1, y = r0_w / log(n), color = "red"), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(0, 2.25), expand = expansion(mult = c(0, 0.01))) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.position = "none"
  ) +
  xlab(expression("W"["1"])) +
  ylab(expression("R" / "log(n)")))

column_label_2 <- wrap_elements(panel = ggpubr::text_grob(
  label = "Pseudo-polar",
  face = "bold",
  size = 15
))
column_label_1 <- wrap_elements(panel = ggpubr::text_grob(
  label = "Euclidean",
  face = "bold",
  size = 15
))
threshold_plots <- (column_label_1 | column_label_2) /
  (data_euc_plot | data_polar_plot) +
  plot_layout(heights = c(0.15, 1)) &
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent")
  )
ggsave("figures/redstone_threshold_plots.pdf",
  plot = threshold_plots,
  dpi = 320,
  bg = "transparent",
  width = 11, height = 5.5
)
knitr::plot_crop("figures/redstone_threshold_plots.pdf")
