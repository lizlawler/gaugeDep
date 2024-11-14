library(evd)
library(mvtnorm)
library(tidyverse)
library(patchwork)

n <- 500
N <- n^2

## Gaussian copula dependence structure ----------
gauss_copula <- function(n = 1000, rho = 0.5) {
  x <- rmvnorm(n, mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  y1 <- qexp(u1)
  y2 <- qexp(u2)
  return(cbind(y1, y2))
}

gauss_gauge <- function(x, y, rho = 0.5) {
  top <- x + y - 2 * rho * sqrt(x * y)
  return(top/(1-rho^2))
}

low_corr_gauss <- gauss_copula(N, 0.1)
mid_corr_gauss <- gauss_copula(N, 0.5)
high_corr_gauss <- gauss_copula(N, 0.9)

gauss_gauge_low <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = gauss_gauge(x, y, 0.1) - 1)
low_corr_gauss_plot <- low_corr_gauss %>% as_tibble() %>% mutate(y1 = y1/log(N), y2 = y2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(gauss_gauge_low$x, gauss_gauge_low$y, z = gauss_gauge_low$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.1,", " expo margins")))
low_corr_gauss_hist <- ggExtra::ggMarginal(low_corr_gauss_plot, type = "histogram")

low_corr_gauss_rw_plot <- low_corr_gauss %>% as_tibble() %>% 
  mutate(r = (y1 + y2), w = y1/r, gw = pmap_dbl(list(x = w, y = 1-w, r = 0.1), function(x, y, r) gauss_gauge(x, y, r))) %>% 
  arrange(w) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw), col = 'red') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.1,", " expo margins")))

gauss_gauge_mid <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = gauss_gauge(x, y, 0.5) - 1)
mid_corr_gauss_plot <- mid_corr_gauss %>% as_tibble() %>% mutate(y1 = y1/log(N), y2 = y2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(gauss_gauge_mid$x, gauss_gauge_mid$y, z = gauss_gauge_mid$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.5,", " expo margins")))
mid_corr_gauss_hist <- ggExtra::ggMarginal(mid_corr_gauss_plot, type = "histogram")

mid_corr_gauss_rw <- mid_corr_gauss %>% as_tibble() %>% 
  mutate(r = (y1 + y2), w = y1/r, gw = pmap_dbl(list(x = w, y = 1-w, r = 0.5), function(x, y, r) gauss_gauge(x, y, r))) %>% 
  arrange(w) 

mid_corr_gauss_rw_plot <- mid_corr_gauss_rw %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw), col = 'red') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.5,", " expo margins")))

gauss_gauge_high <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.9), function(x, y, r) gauss_gauge(x, y, r)-1))
high_corr_gauss_plot <- high_corr_gauss %>% as_tibble() %>% mutate(y1 = y1/log(N), y2 = y2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(gauss_gauge_high$x, gauss_gauge_high$y, z = gauss_gauge_high$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.9,", " expo margins")))
(high_corr_gauss_hist <- ggExtra::ggMarginal(high_corr_gauss_plot, type = "histogram"))

high_corr_gauss_rw_plot <- high_corr_gauss %>% as_tibble() %>% 
  mutate(r = (y1 + y2), w = y1/r, gw = pmap_dbl(list(x = w, y = 1-w, r = 0.9), function(x, y, r) gauss_gauge(x, y, r))) %>% 
  arrange(w) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw), col = 'red') +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.9,", " expo margins")))

all_gauss_plots <- low_corr_gauss_plot + mid_corr_gauss_plot + high_corr_gauss_plot + 
  low_corr_gauss_rw_plot + mid_corr_gauss_rw_plot + high_corr_gauss_rw_plot + plot_layout(nrow = 2)
ggsave("gauss_cops_plots.png", plot = all_gauss_plots, height = 10, width = 15)


## Logistic copula dependence structure -----------
logistic_copula <- function(n = 1000, r = 0.5) {
  x <- rbvevd(n, dep = r, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  y1 <- qexp(u1)
  y2 <- qexp(u2)
  return(cbind(y1, y2))
}

logistic_gauge <- function(x, y, r = 0.5) {
  r_inv <- 1/r
  return(r_inv * pmax(x, y) + (1-r_inv)*pmin(x,y))
}



high_dep_log <- logistic_copula(N, 0.1)
mid_dep_log <- logistic_copula(N, 0.5)
low_dep_log <- logistic_copula(N, 0.9)

logistic_gauge_high <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.1), function(x, y, r) logistic_gauge(x, y, r)-1))
high_dep_log_plot <- high_dep_log %>% as_tibble() %>% mutate(y1 = y1/log(N), y2 = y2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(logistic_gauge_high$x, logistic_gauge_high$y, z = logistic_gauge_high$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Logistic copula, r", "=0.1,", " expo margins")))
(high_dep_log_hist <- ggExtra::ggMarginal(high_dep_log_plot, type = "histogram"))

high_dep_log_rw_plot <- high_dep_log %>% as_tibble() %>% 
  mutate(r = (y1 + y2), 
         w = y1/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, r = 0.1), function(x, y, r) logistic_gauge(x,y,r))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Logistic copula, r", "=0.1,", " expo margins")))

logistic_gauge_mid <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.5), function(x, y, r) logistic_gauge(x, y, r)-1))
mid_dep_log_plot <- mid_dep_log %>% as_tibble() %>% mutate(y1 = y1/log(N), y2 = y2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(logistic_gauge_mid$x, logistic_gauge_mid$y, z = logistic_gauge_mid$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Logistic copula, r", "=0.5,", " expo margins")))
(mid_dep_log_hist <- ggExtra::ggMarginal(mid_dep_log_plot, type = "histogram"))

mid_dep_log_rw_plot <- mid_dep_log %>% as_tibble() %>% 
  mutate(r = (y1 + y2), 
         w = y1/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, r = 0.5), function(x, y, r) logistic_gauge(x,y,r))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Logistic copula, r", "=0.5,", " expo margins")))

logistic_gauge_low <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.9), function(x, y, r) logistic_gauge(x, y, r)-1))
low_dep_log_plot <- low_dep_log %>% as_tibble() %>% mutate(y1 = y1/log(N), y2 = y2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(logistic_gauge_low$x, logistic_gauge_low$y, z = logistic_gauge_low$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Logistic copula, r", "=0.9,", " expo margins")))
(low_dep_log_hist <- ggExtra::ggMarginal(low_dep_log_plot, type = "histogram"))

low_dep_log_rw_plot <- low_dep_log %>% as_tibble() %>% 
  mutate(r = (y1 + y2), 
         w = y1/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, r = 0.9), function(x, y, r) logistic_gauge(x,y,r))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Logistic copula, r", "=0.9,", " expo margins")))

all_log_plots <- low_dep_log_plot + mid_dep_log_plot + high_dep_log_plot + 
  low_dep_log_rw_plot + mid_dep_log_rw_plot + high_dep_log_rw_plot + plot_layout(nrow = 2)
ggsave("logistic_cops_plots.png", plot = all_log_plots, height = 10, width = 15)

## Inverted logistic dependence structure ----------
inv_log_copula <- function(n = 1000, r = 0.5) {
  x <- rbvevd(n, dep = r, model = "log", mar1=c(1,1,1))
  y <- 1/x
  return(y)
}

inv_log_gauge <- function(x, y, r = 0.5) ((x^(1/r) + y^(1/r))^r)

low_dep_invlog <- inv_log_copula(N, 0.9)
mid_dep_invlog <- inv_log_copula(N, 0.5)
high_dep_invlog <- inv_log_copula(N, 0.1)

invlog_gauge_low <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.9), function(x, y, r) inv_log_gauge(x, y, r)-1))
low_dep_invlog_plot <- low_dep_invlog %>% as_tibble() %>% mutate(y1 = V1/log(N), y2 = V2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(invlog_gauge_low$x, invlog_gauge_low$y, z = invlog_gauge_low$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Inv. logistic copula, r", "=0.9,", " expo margins")))
(low_dep_invlog_hist <- ggExtra::ggMarginal(low_dep_invlog_plot, type = "histogram"))

low_dep_invlog_rw_plot <- low_dep_invlog %>% as_tibble() %>% 
  rename(y1 = V1, y2 = V2) %>%
  mutate(r = (y1 + y2), 
         w = y1/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, r = 0.9), function(x, y, r) inv_log_gauge(x,y,r))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +  
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Inv. logistic copula, r", "=0.9,", " expo margins")))

invlog_gauge_mid <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.5), function(x, y, r) inv_log_gauge(x, y, r)-1))
mid_dep_invlog_plot <- mid_dep_invlog %>% as_tibble() %>% mutate(y1 = V1/log(N), y2 = V2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(invlog_gauge_mid$x, invlog_gauge_mid$y, z = invlog_gauge_mid$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Inv. logistic copula, r", "=0.5,", " expo margins")))
(mid_dep_invlog_hist <- ggExtra::ggMarginal(mid_dep_invlog_plot, type = "histogram"))

mid_dep_invlog_rw_plot <- mid_dep_invlog %>% as_tibble() %>% 
  rename(y1 = V1, y2 = V2) %>%
  mutate(r = (y1 + y2), 
         w = y1/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, r = 0.5), function(x, y, r) inv_log_gauge(x,y,r))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +  
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Inv. logistic copula, r", "=0.5,", " expo margins")))

invlog_gauge_high <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.1), function(x, y, r) inv_log_gauge(x, y, r)-1))
high_dep_invlog_plot <- high_dep_invlog %>% as_tibble() %>% mutate(y1 = V1/log(N), y2 = V2/log(N)) %>%
  ggplot(aes(y1, y2)) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(invlog_gauge_high$x, invlog_gauge_high$y, z = invlog_gauge_high$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Inv. logistic copula, r", "=0.1,", " expo margins")))
(high_dep_invlog_hist <- ggExtra::ggMarginal(high_dep_invlog_plot, type = "histogram"))

high_dep_invlog_rw_plot <- high_dep_invlog %>% as_tibble() %>% 
  rename(y1 = V1, y2 = V2) %>%
  mutate(r = (y1 + y2), 
         w = y1/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, r = 0.1), function(x, y, r) inv_log_gauge(x,y,r))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +  
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Inv. logistic copula, r", "=0.1,", " expo margins")))

all_invlog_plots <- low_dep_invlog_plot + mid_dep_invlog_plot + high_dep_invlog_plot + 
  low_dep_invlog_rw_plot + mid_dep_invlog_rw_plot + high_dep_invlog_rw_plot + 
  plot_layout(nrow = 2)
ggsave("invlog_cops_plots.png", plot = all_invlog_plots, height = 10, width = 15)

## Asymmetric logistic dependence structure -------
asym_log_copula <- function(n = 1000, r = 0.5, t1 = 0.5, t2 = 0.5) {
  x <- rbvevd(n, dep = r, model = "alog", asy = c(t1, t2))
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  y1 <- qexp(u1)
  y2 <- qexp(u2)
  return(cbind(y1, y2))
} 

asym_log_gauge <- function(x, y, r = 0.5) {
  r_inv <- 1/r
  return(min((x + y), (r_inv * max(x, y) + (1-r_inv)*min(x,y))))
}

high_dep_asymlog <- asym_log_copula(N, r = 0.1, t1 = 0.9, t2 = 0.5)
mid_dep_asymlog <- asym_log_copula(N, 0.5, t1 = 0.4, t2 = 0.75)
low_dep_asymlog <- asym_log_copula(N, 0.9, t1 = 0.75, t2 = 0.75)

asymlog_gauge_high <- crossing(x = seq(0, 1, length.out = n), y = x) %>% mutate(z = pmap_dbl(list(x = x, y = y, r = 0.1), function(x, y, r) asym_log_gauge(x, y, r)-1))
high_dep_asymlog_plot <- high_dep_asymlog %>% 
  as_tibble() %>% 
  ggplot(aes(y1/log(N), y2/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  # geom_contour(aes(y1, y2, z = 1/g, col = 'red'), breaks = 0) +
  geom_contour(aes(asymlog_gauge_high$x, asymlog_gauge_high$y, z = asymlog_gauge_high$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Asym. logistic copula, r", "=0.1,", " expo margins")))
(high_dep_asymlog_hist <- ggExtra::ggMarginal(high_dep_asymlog_plot, type = "histogram"))

high_dep_asymlog_plot_v2 <- high_dep_asymlog %>% 
  as_tibble() %>% 
  ggplot(aes(y1/log(N), y2/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  # geom_contour(aes(y1, y2, z = 1/g, col = 'red'), breaks = 0) +
  geom_contour(aes(asymlog_gauge_high$x, asymlog_gauge_high$y, z = asymlog_gauge_high$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Asym. logistic copula, r", "=0.1,", " expo margins")))
(high_dep_asymlog_hist_v2 <- ggExtra::ggMarginal(high_dep_asymlog_plot_v2, type = "histogram"))

high_dep_asymlog_plot_rw <- high_dep_asymlog %>% as_tibble() %>% 
  mutate(r = (y1 + y2), 
         w = y1/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, r = 0.1), function(x, y, r) asym_log_gauge(x,y,r))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)") + ggtitle(expression(paste("Asym. logistic copula, r", "=0.1,", " expo margins")))

## Dirichlet dependence structure -------
dirichlet_copula <- function(n = 1000, theta1, theta2) {
  x <- rbvevd(N, alpha = theta1, beta = theta2, model = 'ct')
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y))
}

dirichlet_gauge <- function(x, y, theta1, theta2) {
  return((1 + theta1 + theta2) * pmax(x, y) - (theta1 * x + theta2 * y))
}

w <- seq(0, 1, length.out = 500)
gw <- dirichlet_gauge(w, 1-w, 3, 3)
plot(w/gw, (1-w)/gw)


dirichlet_mid <- dirichlet_copula(N, 7, 2)

dirichlet_mid_gauge <- crossing(x = seq(0, 1, length.out = n), y = x) %>% 
  mutate(z = pmap_dbl(list(x = x, y = y, theta1 = 2, theta2=2), 
                      function(x, y, theta1, theta2) dirichlet_gauge(x, y, theta1, theta2)-1))
# dirichlet_mid_gauge %>% ggplot() + geom_contour(aes(x, y, z = z, col = 'red'), breaks = 0)

dirichlet_mid_plot <- dirichlet_mid %>% 
  as_tibble() %>% 
  ggplot(aes(x/log(N), y/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_contour(aes(dirichlet_mid_gauge$x, dirichlet_mid_gauge$y, z = dirichlet_mid_gauge$z, col = 'red'), breaks = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("x") + ylab("y") + ggtitle("Dirichlet copula")
(dirichlet_mid_plot <- ggExtra::ggMarginal(dirichlet_mid_plot, type = "histogram"))

dirichlet_mid_rw_plot <- dirichlet_mid %>% as_tibble() %>% 
  mutate(r = (x + y), 
         w = x/r, 
         gw = pmap_dbl(list(x = w, y = 1-w, theta1 = 2, theta2 = 2), 
                       function(x, y, theta1, theta2) dirichlet_gauge(x,y,theta1,theta2))) %>%
  ggplot(aes(w, r/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  geom_line(aes(w, 1/gw, col = 'red')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlab("w") + ylab("r/log(N)")

## Husler-Reiss dependence structure
hr_copula <- function(n = 1000, r) {
  x <- rbvevd(N, dep = r, model = 'hr')
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y))
}

hr_gauge <- function(x, y) {
  return(ifelse(x == y, x, Inf))
}
w <- seq(0, 1, length.out = 1000)
gw <- hr_gauge(w, 1-w)
plot(w/gw, (1-w)/gw)

test <- hr_copula(N, 3)
test %>% 
  as_tibble() %>% 
  ggplot(aes(x/log(N), y/log(N))) + 
  geom_point(alpha=0.5, color = "blue") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none")
