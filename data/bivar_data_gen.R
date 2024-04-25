library(evd)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(cmdstanr)
# 
model <- cmdstan_model("stan/bivar_trunc_rectangular.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)

n <- 10000 

## Gaussian dependence structure -----------------------
gauss_copula <- function(n = 1000, rho = 0.5) {
  x <- rmvnorm(n, mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  return(cbind(x, y))
}

low_dep_gauss <- gauss_copula(n, 0.1)
mid_dep_gauss <- gauss_copula(n, 0.5)
high_dep_gauss <- gauss_copula(n, 0.9)

low_dep_gauss_points <- low_dep_gauss %>% as_tibble() %>% 
  mutate(q1 = quantile(low_dep_gauss[,1], 0.95), 
         q2 = quantile(low_dep_gauss[,2], 0.95),
         high = case_when(x >= q1 | y >= q2 ~ 1,
                             .default = 0),
         high = as.factor(high),
         r = x + y,
         w = x/r,
         x_lb = ifelse(w > 0.5, q1, q2 * x/y),
         y_lb = ifelse(w < 0.5, q2, q1 * y / x),
         r0_w = ifelse(w > 0.5, q1/w, q2/(1-w)))

low_dep_gauss_stan <- low_dep_gauss_points %>% filter(high == 1) %>% select(r, w)
stan_low_dep_gauss <- list(R = low_dep_gauss_stan$r, 
                           N = nrow(low_dep_gauss_stan), 
                           W = low_dep_gauss_stan$w, 
                           q = unique(low_dep_gauss_points$qy1),
                           d = 2,
                           iden_mat = matrix(c(1,0,0,1), nrow = 2),
                           rho_mat = matrix(c(0,1,1,0), nrow = 2))
write_stan_json(stan_low_dep_gauss, "data/low_dep_gauss.json")

cols <- c("lightblue", "blue", "red")
low_dep_gauss_plot <- low_dep_gauss_points %>%
  ggplot(aes(x/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.1,", " expo margins")))

low_dep_gauss_plot_rw <- low_dep_gauss_points %>%
  ggplot(aes(y = r/log(n), x = w, color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  ylab("r/log(n)") + xlab("w") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.1,", " expo margins")))

mid_dep_gauss_points <- mid_dep_gauss %>% as_tibble() %>% 
  mutate(qy1 = quantile(mid_dep_gauss[,1], 0.95), qy2 = quantile(mid_dep_gauss[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
mid_dep_gauss_stan <- mid_dep_gauss_points %>% filter(high == 1) %>% select(r, w)
stan_mid_dep_gauss <- list(R = mid_dep_gauss_stan$r, 
                           N = nrow(mid_dep_gauss_stan), 
                           W = mid_dep_gauss_stan$w, 
                           q = unique(mid_dep_gauss_points$qy1),
                           iden_mat = array(c(1,0,0,1), dim = c(2,2)),
                           rho_mat = array(c(0,1,1,0), dim = c(2,2)),
                           d = 2)
write_stan_json(stan_mid_dep_gauss, "data/mid_dep_gauss.json")

mid_dep_gauss_plot <- mid_dep_gauss_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.5,", " expo margins")))

mid_dep_gauss_plot_rw <- mid_dep_gauss_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  xlab("w") + ylab("r/log(n)") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.5,", " expo margins")))

high_dep_gauss_points <- high_dep_gauss %>% as_tibble() %>% 
  mutate(qy1 = quantile(high_dep_gauss[,1], 0.95), qy2 = quantile(high_dep_gauss[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
high_dep_gauss_stan <- high_dep_gauss_points %>% filter(high == 1) %>% select(r, w)
stan_high_dep_gauss <- list(R = high_dep_gauss_stan$r, 
                            N = nrow(high_dep_gauss_stan), 
                            W = high_dep_gauss_stan$w, 
                            q = unique(high_dep_gauss_points$qy1),
                            iden_mat = array(c(1,0,0,1), dim = c(2,2)),
                            rho_mat = array(c(0,1,1,0), dim = c(2,2)),
                            d = 2)
write_stan_json(stan_high_dep_gauss, "data/high_dep_gauss.json")

high_dep_gauss_plot <- high_dep_gauss_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.9,", " expo margins")))

high_dep_gauss_plot_rw <- high_dep_gauss_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  xlab("w") + ylab("r/log(n)") + ggtitle(expression(paste("Gaussian copula, ", rho, "=0.9,", " expo margins")))

all_gauss_data_plots <- low_dep_gauss_plot + mid_dep_gauss_plot + high_dep_gauss_plot + 
  low_dep_gauss_plot_rw + mid_dep_gauss_plot_rw + high_dep_gauss_plot_rw + 
  plot_layout(nrow = 2)
ggsave("gauss_data_plots_threshold.png", plot = all_gauss_data_plots, height = 10, width = 15)


## Logistic dependence structure ------------------------
logistic_copula <- function(n = 1000, r = 0.5) {
  x <- rbvevd(n, dep = r, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  y1 <- qexp(u1)
  y2 <- qexp(u2)
  return(cbind(y1, y2))
}

high_dep_log <- logistic_copula(n, 0.1)
mid_dep_log <- logistic_copula(n, 0.5)
low_dep_log <- logistic_copula(n, 0.9)

high_dep_log_points <- high_dep_log %>% as_tibble() %>% 
  mutate(qy1 = quantile(high_dep_log[,1], 0.95), qy2 = quantile(high_dep_log[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
high_dep_log_stan <- high_dep_log_points %>% filter(high == 1) %>% select(r, w)
stan_high_dep_log <- list(R = high_dep_log_stan$r, 
                          N = nrow(high_dep_log_stan), 
                          W = high_dep_log_stan$w, 
                          q = unique(high_dep_log_points$qy1))
write_stan_json(stan_high_dep_log, "data/high_dep_logistic.json")

high_dep_log_plot <- high_dep_log_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Logistic copula, ", "r", "=0.1,", " expo margins")))

high_dep_log_plot_rw <- high_dep_log_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  xlab("w") + ylab("r/log(n)") + ggtitle(expression(paste("Logistic copula, ", "r", "=0.1,", " expo margins")))

mid_dep_log_points <- mid_dep_log %>% as_tibble() %>% 
  mutate(qy1 = quantile(mid_dep_log[,1], 0.95), qy2 = quantile(mid_dep_log[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
mid_dep_log_stan <- mid_dep_log_points %>% filter(high == 1) %>% select(r, w)
stan_mid_dep_log <- list(R = mid_dep_log_stan$r, 
                         N = nrow(mid_dep_log_stan), 
                         W = mid_dep_log_stan$w, 
                         q = unique(mid_dep_log_points$qy1))
write_stan_json(stan_mid_dep_log, "data/mid_dep_logistic.json")

w <- seq(0, 1, length.out = nrow(mid))

mid_dep_log_plot <- mid_dep_log_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5,size = 2) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent')) +
  scale_x_continuous(limits = c(0,1.2), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1.2), expand = c(0,0)) +
  xlab(expression("X"["1"]/"log(n)")) + ylab(expression("X"["2"]/"log(n)"))

mid_dep_log_plot_rw <- mid_dep_log_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5, size = 2) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        legend.position = "none",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent')) +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,2), expand = c(0,0)) +
  xlab("w") + ylab("r/log(n)")

column_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = 'Pseudo-polar coordinates',
                                                          face = "bold",
                                                          size = 20))
column_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = 'Euclidean coordinates',
                                                          face = "bold",
                                                          size = 20))
threshold_plots <- (column_label_2 | column_label_1) / 
  (mid_dep_log_plot | mid_dep_log_plot_rw) +
  plot_layout(heights = c(0.15, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/prelim_presentation/threshold_plots.pdf",
       plot = threshold_plots,
       dpi = 320,
       bg = 'transparent', 
       width = 11, height = 5.5)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/threshold_plots.pdf")

low_dep_log_points <- low_dep_log %>% as_tibble() %>% 
  mutate(qy1 = quantile(low_dep_log[,1], 0.95), qy2 = quantile(low_dep_log[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
low_dep_log_stan <- low_dep_log_points %>% filter(high == 1) %>% select(r, w)
stan_low_dep_log <- list(R = low_dep_log_stan$r, 
                         N = nrow(low_dep_log_stan), 
                         W = low_dep_log_stan$w, 
                         q = unique(low_dep_log_points$qy1))
write_stan_json(stan_low_dep_log, "data/low_dep_logistic.json")

low_dep_log_plot <- low_dep_log_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Logistic copula, ", "r", "=0.9,", " expo margins")))

low_dep_log_plot_rw <- low_dep_log_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  xlab("w") + ylab("r/log(n)") + ggtitle(expression(paste("Logistic copula, ", "r", "=0.9,", " expo margins")))

all_log_data_plots <- low_dep_log_plot + mid_dep_log_plot + high_dep_log_plot + 
  low_dep_log_plot_rw + mid_dep_log_plot_rw + high_dep_log_plot_rw + 
  plot_layout(nrow = 2)
ggsave("logistic_data_plots_threshold.png", plot = all_log_data_plots, height = 10, width = 15)

## Inverse logistic dependence sturcture ---------------------------------------------
inv_log_copula <- function(n = 1000, r = 0.5) {
  x <- rbvevd(n, dep = r, model = "log", mar1=c(1,1,1))
  y <- 1/x
  return(y)
}

low_dep_invlog <- inv_log_copula(n, 0.9)
mid_dep_invlog <- inv_log_copula(n, 0.5)
high_dep_invlog <- inv_log_copula(n, 0.1)

low_dep_invlog_points <- low_dep_invlog %>% as_tibble() %>% rename(y1 = V1, y2 = V2) %>%
  mutate(qy1 = quantile(low_dep_invlog[,1], 0.95), qy2 = quantile(low_dep_invlog[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
low_dep_invlog_stan <- low_dep_invlog_points %>% filter(high == 1) %>% select(r, w)
stan_low_dep_invlog <- list(R = low_dep_invlog_stan$r, 
                            N = nrow(low_dep_invlog_stan), 
                            W = low_dep_invlog_stan$w, 
                            q = unique(low_dep_invlog_points$qy1))
write_stan_json(stan_low_dep_invlog, "data/low_dep_invlog.json")

low_dep_invlog_plot <- low_dep_invlog_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Inverse logistic copula, ", "r", "=0.9,", " expo margins")))

low_dep_invlog_plot_rw <- low_dep_invlog_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  xlab("w") + ylab("r/log(n)") + ggtitle(expression(paste("Inverse logistic copula, ", "r", "=0.9,", " expo margins")))

mid_dep_invlog_points <- mid_dep_invlog %>% as_tibble() %>% rename(y1 = V1, y2 = V2) %>%
  mutate(qy1 = quantile(mid_dep_invlog[,1], 0.95), qy2 = quantile(mid_dep_invlog[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
mid_dep_invlog_stan <- mid_dep_invlog_points %>% filter(high == 1) %>% select(r, w)
stan_mid_dep_invlog <- list(R = mid_dep_invlog_stan$r, 
                            N = nrow(mid_dep_invlog_stan), 
                            W = mid_dep_invlog_stan$w, 
                            q = unique(mid_dep_invlog_points$qy1))
write_stan_json(stan_mid_dep_invlog, "data/mid_dep_invlog.json")

mid_dep_invlog_plot <- mid_dep_invlog_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Inverse logistic copula, ", "r", "=0.5,", " expo margins")))

mid_dep_invlog_plot_rw <- mid_dep_invlog_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  xlab("w") + ylab("r/log(n)") + ggtitle(expression(paste("Inverse logistic copula, ", "r", "=0.5,", " expo margins")))

high_dep_invlog_points <- high_dep_invlog %>% as_tibble() %>% rename(y1 = V1, y2 = V2) %>%
  mutate(qy1 = quantile(high_dep_invlog[,1], 0.95), qy2 = quantile(high_dep_invlog[,1], 0.95),
         high = case_when(y1 >= qy1 | y2 >= qy2 ~ 1,
                          .default = 0),
         high = as.factor(high),
         r = y1 + y2,
         w = y1/r,
         r0_w = ifelse(w < 0.5, qy1/(1-w), qy2/w),
         y1_lb = ifelse(w < 0.5, qy2, qy1*y2/y1),
         y2_lb = ifelse(w > 0.5, qy1, qy2*y1/y2))
high_dep_invlog_stan <- high_dep_invlog_points %>% filter(high == 1) %>% select(r, w)
stan_high_dep_invlog <- list(R = high_dep_invlog_stan$r, 
                             N = nrow(high_dep_invlog_stan), 
                             W = high_dep_invlog_stan$w, 
                             q = unique(high_dep_invlog_points$qy1))
write_stan_json(stan_high_dep_invlog, "data/high_dep_invlog.json")

high_dep_invlog_plot <- high_dep_invlog_points %>%
  ggplot(aes(y1/log(n), y2/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = y1_lb/log(n), y = y2_lb/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,1) +
  xlab("y1") + ylab("y2") + ggtitle(expression(paste("Inverse logistic copula, ", "r", "=0.1,", " expo margins")))

high_dep_invlog_plot_rw <- high_dep_invlog_points %>%
  ggplot(aes(w, y = r/log(n), color = high)) + 
  geom_point(alpha=0.5) +
  geom_line(aes(x = w, y = r0_w/log(n), color = 'red')) +
  theme_classic() +
  scale_color_manual(values = cols) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14), legend.position = "none") +
  xlim(0,1) + ylim(0,2) +
  xlab("w") + ylab("r/log(n)") + ggtitle(expression(paste("Inverse logistic copula, ", "r", "=0.1,", " expo margins")))

all_invlog_data_plots <- low_dep_invlog_plot + mid_dep_invlog_plot + high_dep_invlog_plot + 
  low_dep_invlog_plot_rw + mid_dep_invlog_plot_rw + high_dep_invlog_plot_rw + 
  plot_layout(nrow = 2)
ggsave("invlog_data_plots_threshold.png", plot = all_invlog_data_plots, height = 10, width = 15)
