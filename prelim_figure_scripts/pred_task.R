# library(evd)
# library(extRemes)
library(mvtnorm)
library(tidyverse)
library(patchwork)

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

e_df <- gauss_copula(1000, 0.9)[[3]]
pred_task_joint <- e_df |> as_tibble() |> ggplot(aes(x = e1, y = e2)) + 
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1))) +
  xlim(0, 10) + ylim(0, 10) +
  annotate("rect", xmin = 8, xmax = 9,   ymin = 8,  ymax = 9, fill = "lightgreen", color = "black", alpha = 1) + 
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

pred_task_x_extreme <- e_df |> as_tibble() |> ggplot(aes(x = e1, y = e2)) + 
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1))) +
  xlim(0, 10) + ylim(0, 10) +
  annotate("rect", xmin = 8, xmax = 9,   ymin = 4,  ymax = 5, fill = "lightgreen", color = "black", alpha = 1) + 
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

pred_task_y_extreme <- e_df |> as_tibble() |> ggplot(aes(x = e1, y = e2)) + 
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1))) +
  xlim(0, 10) + ylim(0, 10) +
  annotate("rect", ymin = 8, ymax = 9, xmin = 4, xmax = 5, fill = "lightgreen", color = "black", alpha = 1) + 
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

# 
# annotate("rect", xmin = 8, xmax = 9,   ymin = 6,  ymax = 7, fill = "lightgreen", color = "black", alpha = 1) + 
#   annotate("rect", xmin = 6, xmax = 7,   ymin = 8,  ymax = 9, fill = "lightgreen", color = "black", alpha = 1) + 

ggsave("~/Desktop/csu/prelim_presentation/pred_task_y.pdf",
       plot = pred_task_y_extreme,
       dpi = 320,
       width = 4,
       height = 4,
       bg = 'transparent')
knitr::plot_crop("~/Desktop/csu/prelim_presentation/pred_task_y.pdf")

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
  xlim(0, 1.5) + ylim(0, 1.5) + 
  geom_point(aes(x = 0.5, y = 0.75), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.5, yend = 0.75)) + 
  geom_point(aes(x = 0.7, y = 0.25), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.7, yend = 0.25)) +
  geom_point(aes(x = 0.15, y = 0.5), size = 2) + geom_segment(aes(x = 0, y = 0, xend = 0.15, yend = 0.5)) +
  xlab("") + ylab("")
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
  xlim(0, 1.5) + ylim(0, 1.5) + 
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


