library(evd)
library(extRemes)
library(mvtnorm)
library(tidyverse)
library(patchwork)

## Gaussian copula dependence structure ----------
gauss_copula <- function(n = 1000, rho = 0.5) {
  x <- rmvnorm(n, mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  f1 <- (-log(u1))^-1
  f2 <- (-log(u2))^-1
  e1 <- qexp(u1)
  e2 <- qexp(u2)
  return(list(original = x, 
              frechet = cbind(f1, f2),
              expo = cbind(e1, e2)))
}

logistic_copula <- function(n = 1000, r = 0.5) {
  x <- rbvevd(n, dep = r, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x1 <- qnorm(u1)
  x2 <- qnorm(u2)
  f1 <- (-log(u1))^-1
  f2 <- (-log(u2))^-1
  e1 <- qexp(u1)
  e2 <- qexp(u2)
  return(list(original = cbind(x1, x2), 
              frechet = cbind(f1, f2),
              expo = cbind(e1, e2)))
}

## n = 100 ----------
## Asymptotic independence (Gaussian)
n1 <- 100
gauss_100 <- gauss_copula(n1, 0.5)
gauss_100_og <- gauss_100$original |> as_tibble()
gauss_100_frechet <- gauss_100$frechet |> as_tibble()
gauss_100_expo <- gauss_100$expo |> as_tibble()

gauss_100_og_plot <- gauss_100_og |> ggplot(aes(x = V1, y = V2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

gauss_100_frechet_plot <- gauss_100_frechet |> ggplot(aes(x = f1/n1, y = f2/n1)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

gauss_100_expo_plot <- gauss_100_expo |> ggplot(aes(x = e1/log(n1), y = e2/log(n1))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

## Asymptotic dependence (Logistic)
logistic_100 <- logistic_copula(n1, 0.5)
logistic_100_og <- logistic_100$original |> as_tibble()
logistic_100_frechet <- logistic_100$frechet |> as_tibble()
logistic_100_expo <- logistic_100$expo |> as_tibble()

logistic_100_og_plot <- logistic_100_og |> ggplot(aes(x = x1, y = x2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

logistic_100_frechet_plot <- logistic_100_frechet |> ggplot(aes(x = f1/n1, y = f2/n1)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

logistic_100_expo_plot <- logistic_100_expo |> ggplot(aes(x = e1/log(n1), y = e2/log(n1))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

row_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = "AI",
                                                       face = "bold",
                                                       size = 20,
                                                       color = 'red'))
row_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = "AD",
                                                       face = "bold",
                                                       size = 20,
                                                       color = 'blue'))
column_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = 'Gaussian margins',
                                                          face = "bold",
                                                          size = 20))
column_sublabel_100 <- wrap_elements(panel = ggpubr::text_grob(label = 'n=100',
                                                          size = 12))
column_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = 'Fréchet margins',
                                                           face = "bold",
                                                           size = 20))
column_label_3 <- wrap_elements(panel = ggpubr::text_grob(label = 'Exponential margins',
                                                          face = "bold",
                                                          size = 20))
blank_space <- wrap_elements(panel = ggpubr::text_grob(label = '',
                                                          size = 5))
n100_plots <- ((blank_space | column_label_1 | column_label_2 | column_label_3) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((blank_space | column_sublabel_100 | column_sublabel_100 | column_sublabel_100) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((row_label_1 | gauss_100_og_plot | gauss_100_frechet_plot | gauss_100_expo_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) /
  ((row_label_2 | logistic_100_og_plot | logistic_100_frechet_plot | logistic_100_expo_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) +
  plot_layout(heights = c(0.15, 0.05, 1, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/prelim_presentation/convergence_n100.png",
       plot = n100_plots,
       dpi = 320,
       bg = "transparent",
       width = 12, height = 7)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/convergence_n100.png")

## n = 1000 -------------
## Asymptotic independence (Gaussian)
n2 <- 1000
gauss_1000 <- gauss_copula(n2, 0.5)
gauss_1000_og <- gauss_1000$original |> as_tibble()
gauss_1000_frechet <- gauss_1000$frechet |> as_tibble()
gauss_1000_expo <- gauss_1000$expo |> as_tibble()

gauss_1000_og_plot <- gauss_1000_og |> ggplot(aes(x = V1, y = V2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

gauss_1000_frechet_plot <- gauss_1000_frechet |> ggplot(aes(x = f1/n2, y = f2/n2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

gauss_1000_expo_plot <- gauss_1000_expo |> ggplot(aes(x = e1/log(n2), y = e2/log(n2))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

## Asymptotic dependence (Logistic)
logistic_1000 <- logistic_copula(n2, 0.5)
logistic_1000_og <- logistic_1000$original |> as_tibble()
logistic_1000_frechet <- logistic_1000$frechet |> as_tibble()
logistic_1000_expo <- logistic_1000$expo |> as_tibble()

logistic_1000_og_plot <- logistic_1000_og |> ggplot(aes(x = x1, y = x2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

logistic_1000_frechet_plot <- logistic_1000_frechet |> ggplot(aes(x = f1/n2, y = f2/n2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

logistic_1000_expo_plot <- logistic_1000_expo |> ggplot(aes(x = e1/log(n2), y = e2/log(n2))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

column_sublabel_1000 <- wrap_elements(panel = ggpubr::text_grob(label = 'n=1000',
                                                               size = 12))
n1000_plots <- ((blank_space | column_label_1 | column_label_2 | column_label_3) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((blank_space | column_sublabel_1000 | column_sublabel_1000 | column_sublabel_1000) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((row_label_1 | gauss_1000_og_plot | gauss_1000_frechet_plot | gauss_1000_expo_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) /
  ((row_label_2 | logistic_1000_og_plot | logistic_1000_frechet_plot | logistic_1000_expo_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) +
  plot_layout(heights = c(0.15, 0.05, 1, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/prelim_presentation/convergence_n1000.png",
       plot= n1000_plots,
       dpi = 320,
       bg = "transparent",
       width = 12, height = 7)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/convergence_n1000.png")

## n = 10000 ---------------------
## Asymptotic independence (Gaussian)
n3 <- 10000
gauss_10000 <- gauss_copula(n3, 0.5)
gauss_10000_og <- gauss_10000$original |> as_tibble()
gauss_10000_frechet <- gauss_10000$frechet |> as_tibble()
gauss_10000_expo <- gauss_10000$expo |> as_tibble()

gauss_10000_og_plot <- gauss_10000_og |> ggplot(aes(x = V1, y = V2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

gauss_10000_frechet_plot <- gauss_10000_frechet |> ggplot(aes(x = f1/n3, y = f2/n3)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

gauss_10000_expo_plot <- gauss_10000_expo |> ggplot(aes(x = e1/log(n3), y = e2/log(n3))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

## Asymptotic dependence (Logistic)
logistic_10000 <- logistic_copula(n3, 0.5)
logistic_10000_og <- logistic_10000$original |> as_tibble()
logistic_10000_frechet <- logistic_10000$frechet |> as_tibble()
logistic_10000_expo <- logistic_10000$expo |> as_tibble()

logistic_10000_og_plot <- logistic_10000_og |> ggplot(aes(x = x1, y = x2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

logistic_10000_frechet_plot <- logistic_10000_frechet |> ggplot(aes(x = f1/n3, y = f2/n3)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

logistic_10000_expo_plot <- logistic_10000_expo |> ggplot(aes(x = e1/log(n3), y = e2/log(n3))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

column_sublabel_10000 <- wrap_elements(panel = ggpubr::text_grob(label = 'n=10000',
                                                                size = 12))
n10000_plots <- ((blank_space | column_label_1 | column_label_2 | column_label_3) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((blank_space | column_sublabel_10000 | column_sublabel_10000 | column_sublabel_10000) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((row_label_1 | gauss_10000_og_plot | gauss_10000_frechet_plot | gauss_10000_expo_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) /
  ((row_label_2 | logistic_10000_og_plot | logistic_10000_frechet_plot | logistic_10000_expo_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) +
  plot_layout(heights = c(0.15, 0.05, 1, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/prelim_presentation/convergence_n10000.png",
       plot = n10000_plots,
       dpi = 320,
       bg = "transparent",
       width = 12, height = 7)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/convergence_n10000.png")

## n = 50000 -------------------------------
## Asymptotic independence (Gaussian)
n4 <- 50000
gauss_50000 <- gauss_copula(n4, 0.5)
gauss_50000_og <- gauss_50000$original |> as_tibble()
gauss_50000_frechet <- gauss_50000$frechet |> as_tibble()
gauss_50000_expo <- gauss_50000$expo |> as_tibble()

gauss_50000_og_plot <- gauss_50000_og |> ggplot(aes(x = V1, y = V2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

gauss_50000_frechet_plot <- gauss_50000_frechet |> ggplot(aes(x = f1/n4, y = f2/n4)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

gauss_50000_expo_plot <- gauss_50000_expo |> ggplot(aes(x = e1/log(n4), y = e2/log(n4))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

gauss_gauge <- function(x, y, rho = 0.5) {
  top <- x + y - 2 * rho * sqrt(x * y)
  return(top/(1-rho^2))
}

w <- seq(0,1, length.out = n4)
gw <- gauss_gauge(w, 1-w, 0.5)
gauss_expo_gauge_plot <- gauss_50000_expo_plot + 
  geom_segment(aes(x = 1, y = 0, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_segment(aes(x = 0, y = 1, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_path(aes(x = w/gw, y = (1-w)/gw), linewidth = 1, color = "blue")

gauss_frechet_rays_plot <- gauss_50000_frechet_plot + 
  geom_abline(intercept = 0.2, slope = -1, linetype = 2, color = 'blue') +
  geom_segment(aes(x = 0, y = 0, xend=gauss_50000_frechet[idx_g[1], ]$f1/n4, yend=gauss_50000_frechet[idx_g[1], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=gauss_50000_frechet[idx_g[3], ]$f1/n4, yend=gauss_50000_frechet[idx_g[3], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=gauss_50000_frechet[idx_g[4], ]$f1/n4, yend=gauss_50000_frechet[idx_g[4], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=gauss_50000_frechet[idx_g[5], ]$f1/n4, yend=gauss_50000_frechet[idx_g[5], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=gauss_50000_frechet[idx_g[6], ]$f1/n4, yend=gauss_50000_frechet[idx_g[6], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=gauss_50000_frechet[idx_g[7], ]$f1/n4, yend=gauss_50000_frechet[idx_g[7], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=gauss_50000_frechet[idx_g[8], ]$f1/n4, yend=gauss_50000_frechet[idx_g[8], ]$f2/n4), linewidth = 1, color = "red")

r_gauss <- gauss_50000_frechet$f1 + gauss_50000_frechet$f2
idx_g <- which(r_gauss/n4 > 0.2)
round(gauss_50000_frechet[idx_g, ]/n4, 2)

## Asymptotic dependence (Logistic)
logistic_50000 <- logistic_copula(n4, 0.5)
logistic_50000_og <- logistic_50000$original |> as_tibble()
logistic_50000_frechet <- logistic_50000$frechet |> as_tibble()
logistic_50000_expo <- logistic_50000$expo |> as_tibble()

logistic_50000_og_plot <- logistic_50000_og |> ggplot(aes(x = x1, y = x2)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(-4, 4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-4, 4), expand = c(0,0)) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

logistic_50000_frechet_plot <- logistic_50000_frechet |> ggplot(aes(x = f1/n4, y = f2/n4)) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  xlab(expression("Z"["1"]/"n")) + ylab(expression("Z"["2"]/"n"))

logistic_50000_expo_plot <- logistic_50000_expo |> ggplot(aes(x = e1/log(n4), y = e2/log(n4))) + geom_point() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_x_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) +
  xlab(expression("Y"["1"]/"log(n)")) + ylab(expression("Y"["2"]/"log(n)"))

rad50000 <- (logistic_50000_frechet$f1 + logistic_50000_frechet$f2)/n4
idx <- which(rad50000 > 0.2)
xpts <- logistic_50000_frechet$f1[idx]
ypts <- logistic_50000_frechet$f2[idx]

logistic_50000_frechet_plot + geom_segment(aes(x = 0, y = 0, xend = xpts[1], yend = ypts[1])) +
  + geom_segment(aes(x = 0, y = 0, xend = xpts[1], yend = ypts[1]))

logistic_50000_expo_plot <- logistic_50000_expo |> ggplot(aes(x = e1/log(n4), y = e2/log(n4))) + geom_point() +
  theme_classic() + 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent')) +
  xlim(0, 1.5) + ylim(0, 1.5) +
  xlab(expression("X"["1"]/"log(n)")) + ylab(expression("X"["2"]/"log(n)"))

logistic_gauge <- function(x, y, r = 0.5) {
  r_inv <- 1/r
  return(r_inv * pmax(x, y) + (1-r_inv)*pmin(x,y))
}
gw_log <- logistic_gauge(w, 1-w, 0.5)

logistic_expo_gauge_plot <- logistic_50000_expo_plot + 
  geom_segment(aes(x = 1, y = 0, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_segment(aes(x = 0, y = 1, xend=1, yend=1), linetype = 2, color = "lightgrey") +
  geom_path(aes(x = w/gw_log, y = (1-w)/gw_log), linewidth = 1, color = "blue")
  

r_logistic <- logistic_50000_frechet$f1 + logistic_50000_frechet$f2
idx <- which(r_logistic/n4 > 0.2)
idx <- idx[c(1:4, 6:8, 10)]
round(logistic_50000_frechet[idx, ]/n4, 2)

logistic_frechet_rays_plot <- logistic_50000_frechet_plot + 
  geom_abline(intercept = 0.2, slope = -1, linetype = 2, color = 'blue') +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[1], ]$f1/n4, yend=logistic_50000_frechet[idx[1], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[2], ]$f1/n4, yend=logistic_50000_frechet[idx[2], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[3], ]$f1/n4, yend=logistic_50000_frechet[idx[3], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[4], ]$f1/n4, yend=logistic_50000_frechet[idx[4], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[5], ]$f1/n4, yend=logistic_50000_frechet[idx[5], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[6], ]$f1/n4, yend=logistic_50000_frechet[idx[6], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[7], ]$f1/n4, yend=logistic_50000_frechet[idx[7], ]$f2/n4), linewidth = 1, color = "red") +
  geom_segment(aes(x = 0, y = 0, xend=logistic_50000_frechet[idx[8], ]$f1/n4, yend=logistic_50000_frechet[idx[8], ]$f2/n4), linewidth = 1, color = "red")

column_sublabel_50000 <- wrap_elements(panel = ggpubr::text_grob(label = 'n=50000',
                                                                 size = 12))
n50000_plots <- ((blank_space | column_label_1 | column_label_2 | column_label_3) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((blank_space | column_sublabel_50000 | column_sublabel_50000 | column_sublabel_50000) + plot_layout(widths = c(0.15, .95, .95, .95))) / 
  ((row_label_1 | gauss_50000_og_plot | gauss_frechet_rays_plot | gauss_expo_gauge_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) /
  ((row_label_2 | logistic_50000_og_plot | logistic_frechet_rays_plot | logistic_expo_gauge_plot) + plot_layout(widths = c(0.15, .95, .95, .95))) +
  plot_layout(heights = c(0.15, 0.05, 1, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/prelim_presentation/convergence_n50000_rays_gauge.png",
       plot = n50000_plots,
       dpi = 320,
       bg = "transparent",
       width = 12, height =7)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/convergence_n50000_rays_gauge.png")

chi_gauss <- function(u, rho = 0.5) {
  corr <- diag(2)
  corr[lower.tri(corr)] <- rho
  corr[upper.tri(corr)] <- rho
  num <- log(mvtnorm::pmvnorm(lower = -Inf, upper = c(qnorm(u), qnorm(u)), corr = corr)[1])
  denom <- log(u)
  return(2 - num/denom)
}

chi_gauss_list <- lapply(seq(0.85, 0.999999, length.out = 500), function(x) chi_gauss(x, .5))
chi_gauss_df <- chi_gauss_list |> unlist()
plot(x = seq(0.85, 0.999999, length.out = 500), y = chi_gauss_df, type = "l", xlim = c(0.9, 1), ylim = c(0,1))

chi_logistic <- function(u, r = 0.5) {
  gumbel_quant <- qgev(u, loc = 0, scale = 1, shape = 0)
  num <- log(pbvevd(q = c(gumbel_quant, gumbel_quant), dep = r, model = "log"))
  denom <- log(u)
  return(2 - num/denom)
}

chi_and_ratio <- function(u, dep) {
  vol_comp <- 1 - dep
  chi <- chi_logistic(u, dep)
  return(list(chi = chi, ratio = chi/vol_comp))
}
chi_and_ratio(0.95, 0.95)

chi_logistic(0.95, 0.1) / 0.9

pts <- rbvevd(1000, dep = 0.1, model = "log")
pts1_cdf <- pgev()