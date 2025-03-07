library(evd)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(qs)

cols <- c("lightblue", "blue", "red")
data <- qread("data/redstone_expo.qs")
head(data$cloud_tib)
cloud_tib <- data$cloud_tib
n <- nrow(cloud_tib)

(data_euc_plot<- cloud_tib |>
  ggplot(aes(x/log(n), y/log(n), color = high)) + 
  geom_point(alpha=0.8) +
  geom_line(aes(x = x_lb/log(n), y = y_lb/log(n), color = 'red'), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  xlab(expression("ERC"/"log(n)")) + ylab(expression("FWI"/"log(n)")))

(data_polar_plot <- cloud_tib |>
  ggplot(aes(w1, r/log(n), color = high)) + 
  geom_point(alpha=0.8) +
  geom_line(aes(x = w1, y = r0_w/log(n), color = 'red'), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  theme(panel.background = element_rect(fill = "transparent", color = "transparent"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  xlab(expression("W"["1"])) + ylab(expression("R"/"log(n)")))

column_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = 'Pseudo-polar',
                                                          face = "bold",
                                                          size = 15))
column_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = 'Euclidean',
                                                          face = "bold",
                                                          size = 15))
threshold_plots <- (column_label_1 | column_label_2) / 
  (data_euc_plot | data_polar_plot) +
  plot_layout(heights = c(0.15, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/research/posters-presentations/NSA_talk/redstone_threshold_plots.pdf",
       plot = threshold_plots,
       dpi = 320,
       bg = 'transparent', 
       width = 11, height = 5.5)
knitr::plot_crop("~/Desktop/research/posters-presentations/NSA_talk/redstone_threshold_plots.pdf")
