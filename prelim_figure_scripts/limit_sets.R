# library(evd)
# library(extRemes)
library(mvtnorm)
library(tidyverse)
library(patchwork)


gauss_gauge <- function(x, y, rho = 0.5) {
  top <- x + y - 2 * rho * sqrt(x * y)
  return(top/(1-rho^2))
}

w <- seq(0, 1, length.out = 1000)
gw <- gauss_gauge(w, 1-w, 0.5)
df <- cbind(w, 1-w, gw) |> as_tibble() |> rename(v = V2)
limit_set_plot <- df |> ggplot() + 
  geom_segment(aes(x = 1, y = 0, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_segment(aes(x = 0, y = 1, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_path(aes(x = w/gw, y =v/gw), linewidth = 1, color = "blue") +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent')) +
  geom_point(aes(x = 0.5, y = 0.75), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.5, yend = 0.75)) + 
  geom_point(aes(x = 0.7, y = 0.25), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.7, yend = 0.25)) +
  geom_point(aes(x = 0.15, y = 0.5), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.15, yend = 0.5)) +
  xlab("") + ylab("") +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0))
ggsave("~/Desktop/csu/prelim_presentation/limit_set_plot.pdf",
       plot = limit_set_plot,
       dpi = 320,
       bg = "transparent",
       width =3, height = 3)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/limit_set_plot.pdf")


limit_set_geometry <- df |> ggplot() + 
  geom_segment(aes(x = 1, y = 0, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_segment(aes(x = 0, y = 1, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_path(aes(x = w/gw, y =v/gw), linewidth = 1, color = "blue") +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent')) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  geom_segment(aes(x = 0, y = 0, xend = 1.25, yend = 1.25),
               arrow = arrow(length = unit(0.5, "cm"))) +
  annotate("rect", xmin = 0.75, xmax = Inf,   ymin = 0.75,  ymax = Inf, fill = "salmon", alpha = 0.3) +
  geom_point(aes(x = 0.75, y = 0.75), shape = 8, size = 4, color = "red") +
  geom_point(aes(x = 0.75, y = 0.75), shape = 19, size = 2.5, color = "red") +
  xlab("") + ylab("")
  
# + geom_segment(aes(x = 0, y = 0, xend = 0.5, yend = 0.75)) + 
#   geom_point(aes(x = 0.7, y = 0.25), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.7, yend = 0.25)) +
#   geom_point(aes(x = 0.15, y = 0.5), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.15, yend = 0.5)) +
#   xlab("") + ylab("")
ggsave("~/Desktop/csu/prelim_presentation/limit_set_geometry.pdf",
       plot = limit_set_geometry,
       dpi = 320,
       bg = "transparent",
       width =3, height = 3)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/limit_set_geometry.pdf")


