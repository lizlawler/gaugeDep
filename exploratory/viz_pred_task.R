# =============================================================================
# viz_pred_task.R
#
# Creates the prediction task illustration figure showing the three prediction
# boxes (b1, b2, b3) overlaid on simulated bivariate data under the Gaussian
# and logistic dependence structures. Used to explain the prediction task setup.
#
# Inputs:    (generates data internally)
# Outputs:   figures/gauss_logistic_pred_task.png, figures/husler_reiss_pred_task.png
# =============================================================================

library(patchwork)
library(grafify)
library(tidyverse)
library(evd)
library(mvtnorm)
library(tidyverse)
library(gaugeDependence)

gauss <- function(N = 5000, dep = 0.5) {
  x <- rmvnorm(N, mean = c(0,0), sigma = matrix(c(1, dep, dep, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x,y))
}

logistic <- function(N = 5000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x,y))
}

gauss_low <- gauss(5000, 0.1) |> as_tibble() |> mutate(dep_type = "gauss_low")
gauss_mid <- gauss(5000, 0.5) |> as_tibble() |> mutate(dep_type = "gauss_mid")
gauss_high <- gauss(5000, 0.9) |> as_tibble() |> mutate(dep_type = "gauss_high")
logistic_low <- logistic(5000, 0.9) |> as_tibble() |> mutate(dep_type = "logistic_low")
logistic_mid <- logistic(5000, 0.5) |> as_tibble() |> mutate(dep_type = "logistic_mid")
logistic_high <- logistic(5000, 0.1) |> as_tibble() |> mutate(dep_type = "logistic_high")

all_samps <- rbind(gauss_low, gauss_mid, gauss_high, logistic_low, logistic_mid, logistic_high) |>
  mutate(dep_type = factor(dep_type, levels = c("gauss_low", "gauss_mid", "gauss_high", "logistic_low", "logistic_mid", "logistic_high")))

b1_col <- get_graf_colours("contrast_blue")
b2_col <- get_graf_colours("contrast_red")
b3_col <- get_graf_colours("contrast_yellow")

all_samps |> as_tibble() |> ggplot(aes(x = x, y = y)) + 
  # geom_point(size = 0.5) +
  geom_point(alpha = 0.75) +
  facet_wrap(. ~ dep_type, axes = "all") +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.5)),
        strip.text.x = element_blank(),
        strip.background.x = element_rect(fill='transparent', color = 'black'),
        panel.spacing.x = unit(0.85, "cm", data = NULL),
        panel.spacing.y = unit(0.85, "cm", data = NULL)) +
  scale_x_continuous(limits = c(NA, 13), expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(limits = c(NA, 13), expand = expansion(mult = c(0,0.01))) +
  annotate("rect", xmin = 10, xmax = 12, ymin = 10,  ymax = 12, fill = b1_col, color = "black", alpha = 0.75) +
  annotate("rect", xmin = 10, xmax = 12, ymin = 6,  ymax = 8, fill = b2_col, color = "black", alpha = 0.75) +
  annotate("rect", xmin = 10, xmax = 12, ymin = 2,  ymax = 4, fill = b3_col, color = "black", alpha = 0.75) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

ggsave("figures/gauss_logistic_pred_task.png",
       dpi = 320,
       width = 12,
       height = 8,
       bg = 'transparent')
knitr::plot_crop("figures/gauss_logistic_pred_task.png")

# dep_levels <- list(c(0.1, "low"), c(1, "mid"), c(3, "high"))
husler_reiss <- function(N = 5000, dep = 1) {
  x <- rbvevd(N, dep = dep, model = "hr")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y))
}

hr_low <- husler_reiss(dep = 0.1) |> as_tibble() |> mutate(level = "low")
hr_mid <- husler_reiss(dep = 1) |> as_tibble() |> mutate(level = "mid")
hr_high <- husler_reiss(dep = 3) |> as_tibble() |> mutate(level = "high")
all_samps <- rbind(hr_low, hr_mid, hr_high) |> mutate(level = factor(level, levels = c("low", "mid", "high")))

all_samps |> as_tibble() |> ggplot(aes(x = x, y = y)) + 
  # geom_point(size = 0.5) +
  geom_point(alpha = 0.75) +
  facet_wrap(. ~ level, axes = "all") +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.5)),
        strip.text.x = element_blank(),
        strip.background.x = element_rect(fill='transparent', color = 'black'),
        panel.spacing.x = unit(0.85, "cm", data = NULL),
        panel.spacing.y = unit(0.85, "cm", data = NULL)) +
  scale_x_continuous(limits = c(NA, 13), expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(limits = c(NA, 13), expand = expansion(mult = c(0,0.01))) +
  annotate("rect", xmin = 10, xmax = 12, ymin = 10,  ymax = 12, fill = b1_col, color = "black", alpha = 0.75) +
  annotate("rect", xmin = 10, xmax = 12, ymin = 6,  ymax = 8, fill = b2_col, color = "black", alpha = 0.75) +
  annotate("rect", xmin = 10, xmax = 12, ymin = 2,  ymax = 4, fill = b3_col, color = "black", alpha = 0.75) +
  xlab(expression("X"["1"])) + ylab(expression("X"["2"]))

ggsave("figures/husler_reiss_pred_task.png",
       dpi = 320,
       width = 12,
       height = 4,
       bg = 'transparent')
knitr::plot_crop("figures/husler_reiss_pred_task.png")
