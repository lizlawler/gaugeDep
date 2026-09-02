# =============================================================================
# Visualises the marginal gauge threshold r0(W) as a function of angle W
# for a representative simulated dataset. Shows both Euclidean and polar
# coordinate representations with above-threshold exceedances highlighted.
#
# Inputs:    data/{dep_type}/{dep_level}_{i}.json
# Outputs:   figures/marg_threshold.png
# =============================================================================

library(evd)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(qs)

cols <- c("lightblue", "red", "blue")
data <- RcppSimdJson::fload("data/gauss/mid_100.json")
w <- data$W
r <- data$R
x <- w * r
y <- (1-w)*r
cloud_tib <- tibble(x = x, y=y, w = w, r = r, r0w = data$r0_w) |> 
  mutate(high = r > r0w, 
         x_lb = w*r0w,
         y_lb = (1-w)*r0w)
n <- nrow(cloud_tib)

# ylab(expression("log(T"["est"]~")")) +
#   xlab(expression("log(T"["true"]~")")) +

(data_euc_plot<- cloud_tib |>
    ggplot(aes(x/log(n), y/log(n), color = high)) + 
    geom_point(alpha=0.8) +
    geom_line(aes(x = x_lb/log(n), y = y_lb/log(n), color = 'red'), linewidth = 0.75) +
    theme_classic() +
    scale_color_manual(values = cols) +
    scale_x_continuous(limits = c(NA,1.1), , expand = expansion(mult = c(0,0.01))) +
    scale_y_continuous(limits = c(NA,1.1), expand = expansion(mult = c(0,0.01))) +
    theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          axis.text = element_text(size = rel(1.3)),
          axis.title = element_text(size = rel(1.3)),
          legend.position = "none") +
    xlab(expression("Y"["1"]~"/log(n)")) +
    ylab(expression("Y"["2"]~"/log(n)")))

(data_polar_plot <- cloud_tib |>
    ggplot(aes(w, r/log(n), color = high)) + 
    geom_point(alpha=0.8) +
    geom_line(aes(x = w, y = r0w/log(n), color = 'red'), linewidth = 0.75) +
    theme_classic() +
    scale_color_manual(values = cols) +
    scale_x_continuous(limits = c(NA, 1), expand = expansion(mult = c(0,0.02))) +
    scale_y_continuous(limits = c(NA, 2), expand = expansion(mult = c(0,0.01))) +
    theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
          plot.background = element_rect(fill = "transparent", color = "transparent"),
          axis.text = element_text(size = rel(1.3)),
          axis.title = element_text(size = rel(1.3)),
          legend.position = "none") +
    xlab(expression("W"["1"])) + ylab(expression("R"/"log(n)")))

column_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = 'Pseudo-polar',
                                                          face = "bold",
                                                          size = 15))
column_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = 'Euclidean',
                                                          face = "bold",
                                                          size = 15))
(column_label_1 | column_label_2) / 
  (data_euc_plot | data_polar_plot) +
  plot_layout(heights = c(0.15, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("figures/marg_threshold.png",
       dpi = 320,
       bg = 'transparent', 
       width = 11, height = 5.5)
knitr::plot_crop("figures/marg_threshold.png")
