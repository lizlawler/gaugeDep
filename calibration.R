library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)

quantile_df <- function(x, probs = c(0.25, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = c('lower', 'upper')
  ) |> pivot_wider(names_from = quant, values_from = val)
}

create_fit_df <- function(gauge, dep_level, lhood_type, thresh_type, data_num) {
  csvfiles <- paste0("stan/csv_fits/calibrate/", gauge, "/",
                     list.files(path = paste0("stan/csv_fits/calibrate/", gauge, "/"), 
                                pattern = paste0(dep_level, "_", data_num, "_", lhood_type, "_", thresh_type, "_\\d{1}.csv")))
  fit <- as_cmdstan_fit(csvfiles)
  if (gauge != "dirichlet") {
    return(fit |> as_draws_df() |>
             select(-c(".iteration", ".chain")) |>
             rename(draw = ".draw") |> 
             select(dep))
  } else {
    return(fit |> as_draws_df() |>
             select(-c(".iteration", ".chain")) |>
             rename(draw = ".draw") |> 
             select(theta1, theta2))
  }
}

create_tib_ci <- function(fit_df, ci_level, truth) {
  alpha <- (1-ci_level)/2
  probs_vec <- c(alpha, 1 - alpha)
  if (ncol(fit_df) == 1) {
    temp <- quantile_df(fit_df$dep, probs = probs_vec)
    return(temp |> mutate(level = ci_level, 
                          truth = truth, 
                          coverage = ifelse(lower <= truth & upper >= truth, 1, 0)))
  } else {
    temp <- apply(fit_df, 2, quantile_df, probs = probs_vec) |> 
      bind_rows() |> 
      mutate(param = c("theta1", "theta2"), .before = "lower")
    return(temp |> mutate(level = rep(ci_level, 2), 
                          truth = c(truth, 2), 
                          coverage = ifelse(lower <= truth & upper >= truth, 1, 0)))
  }
}
# 
# # for dirichlet
# dep_levels <- list(c(3, "high"), c(1, "mid"), c(0.5, "low"))

coverage_list <- function(gauge, dep_level, lhood_type, thresh_type, data_num, truth) {
  fit_df <- create_fit_df(gauge, dep_level, lhood_type, thresh_type, data_num)
  levels_vec <- seq(0.05, 0.95, by = 0.05)
  cov_list <- sapply(levels_vec, function(x) create_tib_ci(fit_df, x, truth), simplify = FALSE) |> 
    bind_rows() |>
    mutate(dataset = data_num)
  return(list(fit_df = fit_df, cov_list = cov_list))
}

extract_coverage <- function(gauge, lhood_type, thresh_type) {
  temp_list <- vector("list", 3)
  names(temp_list) <- c("high", "mid", "low")
  if (gauge == "gauss") {
    true_vals <- c(0.9, 0.5, 0.1)
    for(i in seq_along(true_vals)) {
      temp_list[[i]] <- sapply(1:100, function(y) coverage_list("gauss", names(temp_list)[i], lhood_type, thresh_type, y, true_vals[i]))
    }
  } else if (gauge == "dirichlet") {
    true_vals <- c(3, 1, 0.5)  
    for(i in seq_along(true_vals)) {
      temp_list[[i]] <- sapply(1:100, function(y) coverage_list("dirichlet", names(temp_list)[i], lhood_type, thresh_type, y, true_vals[i]))
    }
  } else {
    true_vals <- c(0.1, 0.5, 0.9)
    for(i in seq_along(true_vals)) {
      temp_list[[i]] <- sapply(1:100, function(y) coverage_list(gauge, names(temp_list)[i], lhood_type, thresh_type, y, true_vals[i]))
    }  
  }
  cov_plot_df <- lapply(temp_list, function(x) x[2,] |> bind_rows()) |>
    bind_rows() |>
    group_by(level, truth) |>
    summarize(p_hat = mean(coverage),
              sd = sd(coverage)) |>
    ungroup() |>
    mutate(se = sd/10,
           lb = p_hat - qnorm(0.975) * se,
           ub = p_hat + qnorm(0.975) * se,
           truth = as.factor(truth))
  med_mean_df <- lapply(temp_list, function(x) lapply(x[1,], function(y)  c("median" = median(y$dep), "mean" = mean(y$dep)))) |> bind_rows() |> mutate(truth = rep(true_vals, times = rep(100,3)))
  return(list(cov_plot_df = cov_plot_df,
              med_mean_df = med_mean_df))
}


## Extract coverage from fits -----------
gauss_coverage_marg <- extract_coverage("gauss", "trunc", "marg")
logistic_coverage_marg <- extract_coverage("logistic", "trunc", "marg")

gauss_coverage_ctau <- extract_coverage("gauss", "trunc", "ctau")
logistic_coverage_ctau <- extract_coverage("logistic", "trunc", "ctau")

gauss_coverage_cens <- extract_coverage("gauss", "cens", "marg")


## Create boxplots and coverage plots ----------
plot_coverage <- function(cov_tibble, true_dep) {
  cov_tibble |> filter(truth == true_dep) |> ggplot(aes(x = level, y = p_hat)) + 
    geom_point() +
    geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
    geom_abline() +
    theme_classic() +
    ylim(-0.01,1.01) +
    xlim(-0.01,1.01) +
    xlab("Nominal Rate") +
    ylab("Empirical Rate") +
    theme(panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.text = element_text(size = rel(1.2)),
          axis.title = element_text(size = rel(1.2))) 
}

plot_boxplot <- function(cov_tibble, true_dep) {
  cov_tibble |> filter(truth == true_dep) |> 
    ggplot(aes(y = median, x = 0)) + 
    geom_boxplot(fill = "grey", width = 0.5) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      axis.text = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.2))) +
    geom_hline(yintercept = true_dep, col = "blue", linetype = 2, size = 1.2) +
    ylim(0, 1) +
    xlim(-0.5, 0.5)
}

create_save_plots <- function(coverage_list) {
  
}

all_gauss_marg_plots <- (plot_coverage(gauss_coverage_marg[[1]], 0.1) | plot_coverage(gauss_coverage_marg[[1]], 0.5) | plot_coverage(gauss_coverage_marg[[1]], 0.9)) / 
  (plot_boxplot(gauss_coverage_marg[[2]], 0.1) | plot_boxplot(gauss_coverage_marg[[2]], 0.5) | plot_boxplot(gauss_coverage_marg[[2]], 0.9))

ggsave("./figures/gauss_marg_calibration.pdf", 
       all_gauss_marg_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

all_gauss_ctau_plots <- (plot_coverage(gauss_coverage_ctau[[1]], 0.1) | plot_coverage(gauss_coverage_ctau[[1]], 0.5) | plot_coverage(gauss_coverage_ctau[[1]], 0.9)) / 
  (plot_boxplot(gauss_coverage_ctau[[2]], 0.1) | plot_boxplot(gauss_coverage_ctau[[2]], 0.5) | plot_boxplot(gauss_coverage_ctau[[2]], 0.9))

ggsave("./figures/gauss_ctau_calibration.pdf", 
       all_gauss_ctau_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

all_logistic_marg_plots <- (plot_coverage(logistic_coverage_marg[[1]], 0.1) | plot_coverage(logistic_coverage_marg[[1]], 0.5) | plot_coverage(logistic_coverage_marg[[1]], 0.9)) / 
  (plot_boxplot(logistic_coverage_marg[[2]], 0.1) | plot_boxplot(logistic_coverage_marg[[2]], 0.5) | plot_boxplot(logistic_coverage_marg[[2]], 0.9))

ggsave("./figures/logistic_marg_calibration.pdf", 
       all_logistic_marg_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

all_logistic_ctau_plots <- (plot_coverage(logistic_coverage_ctau[[1]], 0.1) | plot_coverage(logistic_coverage_ctau[[1]], 0.5) | plot_coverage(logistic_coverage_ctau[[1]], 0.9)) / 
  (plot_boxplot(logistic_coverage_ctau[[2]], 0.1) | plot_boxplot(logistic_coverage_ctau[[2]], 0.5) | plot_boxplot(logistic_coverage_ctau[[2]], 0.9))

ggsave("./figures/logistic_ctau_calibration.pdf", 
       all_logistic_ctau_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)





library(patchwork)
column_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = 'Point estimate',
                                                          face = "bold",
                                                          size = 18))
row_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = 'Gaussian',
                                                       face = "bold",
                                                       size = 18))
column_label_3 <- wrap_elements(panel = ggpubr::text_grob(label = 'Calibration',
                                                          face = "bold",
                                                          size = 18))
column_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = '        ',
                                                          size = 5))
row_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = 'Logistic',
                                                       face = "bold",
                                                       size = 18))
calib_plots <- ((column_label_1 | column_label_2 |  column_label_3) + plot_layout(widths = c(0.4, 1, 1))) / 
  ((row_label_1 | (gauss_high_med_boxplot + theme(plot.margin = unit(c(0,50,0,0), "pt")))  |  gauss_calibration_high) + plot_layout(widths = c(0.4, 1, 1))) / 
  ((row_label_2 | (logistic_mid_med_boxplot + theme(plot.margin = unit(c(0,50,0,0), "pt"))) | logistic_calibration_mid) + plot_layout(widths = c(0.4, 1, 1)))  +
  plot_layout(heights = c(.15, 1, 1), width = c(0.01, 1, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/prelim_presentation/calibration_plots.pdf",
       plot = calib_plots,
       dpi = 320,
       bg = 'transparent',
       width = 12,
       height = 8)
knitr::plot_crop("~/Desktop/csu/prelim_presentation/calibration_plots.pdf")

 ## Inverse logistic dependence -----
inv_log_coverage <- inv_log_high_list[2,] |> 
  bind_rows() |> 
  rbind(inv_log_mid_list[2,] |> 
          bind_rows()) |>
  rbind(inv_log_low_list[2,] |> 
          bind_rows()) |>
  group_by(level, truth) |>
  summarize(p_hat = mean(coverage),
            sd = sd(coverage)) |>
  ungroup() |>
  mutate(se = sd/10,
         lb = p_hat - qnorm(0.975) * se,
         ub = p_hat + qnorm(0.975) * se,
         truth = as.factor(truth))

inv_log_calibration_low <- inv_log_coverage |> filter(truth == 0.9) |> ggplot(aes(x = level, y = p_hat)) + 
  geom_point() +
  geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
  geom_abline() +
  theme_classic() +
  ylim(-0.001,1.01) +
  xlim(-0.001,1.01) +
  xlab("Nominal Rate") +
  ylab("Empirical Rate") +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) 
ggsave("./figures/inv_log_calibration_low.pdf", 
       inv_log_calibration_low,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/inv_log_calibration_low.pdf")

inv_log_calibration_mid <- inv_log_coverage |> filter(truth == 0.5) |> ggplot(aes(x = level, y = p_hat)) + 
  geom_point() +
  geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
  geom_abline() +
  theme_classic() +
  ylim(-0.01,1.01) +
  xlim(-0.01,1.01) +
  xlab("Nominal Rate") +
  ylab("Empirical Rate") +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) 
ggsave("./figures/inv_log_calibration_mid.pdf", 
       inv_log_calibration_mid,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/inv_log_calibration_mid.pdf")

inv_log_calibration_high <- inv_log_coverage |> filter(truth == 0.1) |> ggplot(aes(x = level, y = p_hat)) + 
  geom_point() +
  geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
  geom_abline() +
  theme_classic() +
  ylim(-0.01,1.01) +
  xlim(-0.01,1.01) +
  xlab("Nominal Rate") +
  ylab("Empirical Rate") +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) 
ggsave("./figures/inv_log_calibration_high.pdf", 
       inv_log_calibration_high,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/inv_log_calibration_high.pdf")

inv_log_low_med <- lapply(inv_log_low_list[1,], function(x) median(x$dep)) |> 
  unlist() |> 
  as_tibble()
inv_log_low_med_boxplot <- inv_log_low_med |> ggplot(aes(y = value, x = 0)) + 
  geom_boxplot(fill = "grey", width = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)) +
  geom_hline(yintercept = 0.9, col = "blue", linetype = 2) +
  ylim(0, 1) +
  xlim(-0.5, 0.5)
ggsave("./figures/inv_log_low_med_boxplot.pdf", 
       inv_log_low_med_boxplot,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/inv_log_low_med_boxplot.pdf")

inv_log_mid_med <- lapply(inv_log_mid_list[1,], function(x) median(x$dep)) |> 
  unlist() |> 
  as_tibble()
inv_log_mid_med_boxplot <- inv_log_mid_med |> ggplot(aes(y = value, x = 0)) + 
  geom_boxplot(fill = "grey", width = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)) +
  geom_hline(yintercept = 0.5, col = "blue", linetype = 2) +
  ylim(0, 1) +
  xlim(-0.5, 0.5)
ggsave("./figures/inv_log_mid_med_boxplot.pdf", 
       inv_log_mid_med_boxplot,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/inv_log_mid_med_boxplot.pdf")

inv_log_high_med <- lapply(inv_log_high_list[1,], function(x) median(x$dep)) |> 
  unlist() |> 
  as_tibble()
inv_log_high_med_boxplot <- inv_log_high_med |> ggplot(aes(y = value, x = 0)) + 
  geom_boxplot(fill = "grey", width = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)) +
  geom_hline(yintercept = 0.1, col = "blue", linetype = 2) +
  ylim(0, 1) +
  xlim(-0.5, 0.5)
ggsave("./figures/inv_log_high_med_boxplot.pdf", 
       inv_log_high_med_boxplot,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/inv_log_high_med_boxplot.pdf")


## Asymmetric logistic dependence -----
asym_log_coverage <- asym_log_high_list[2,] |> 
  bind_rows() |> 
  rbind(asym_log_mid_list[2,] |> 
          bind_rows()) |>
  rbind(asym_log_low_list[2,] |> 
          bind_rows()) |>
  group_by(level, truth) |>
  summarize(p_hat = mean(coverage),
            sd = sd(coverage)) |>
  ungroup() |>
  mutate(se = sd/10,
         lb = p_hat - qnorm(0.975) * se,
         ub = p_hat + qnorm(0.975) * se,
         truth = as.factor(truth))

asym_log_calibration_low <- asym_log_coverage |> filter(truth == 0.9) |> ggplot(aes(x = level, y = p_hat)) + 
  geom_point() +
  geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
  geom_abline() +
  theme_classic() +
  ylim(-0.1,1.01) +
  xlim(-0.1,1.01) +
  xlab("Nominal Rate") +
  ylab("Empirical Rate") +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) 
ggsave("./figures/asym_log_calibration_low.pdf", 
       asym_log_calibration_low,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/asym_log_calibration_low.pdf")

asym_log_calibration_mid <- asym_log_coverage |> filter(truth == 0.5) |> ggplot(aes(x = level, y = p_hat)) + 
  geom_point() +
  geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
  geom_abline() +
  theme_classic() +
  ylim(-0.01,1.01) +
  xlim(-0.01,1.01) +
  xlab("Nominal Rate") +
  ylab("Empirical Rate") +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) 
ggsave("./figures/asym_log_calibration_mid.pdf", 
       asym_log_calibration_mid,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/asym_log_calibration_mid.pdf")

asym_log_calibration_high <- asym_log_coverage |> filter(truth == 0.1) |> ggplot(aes(x = level, y = p_hat)) + 
  geom_point() +
  geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
  geom_abline() +
  theme_classic() +
  ylim(-0.01,1.01) +
  xlim(-0.01,1.01) +
  xlab("Nominal Rate") +
  ylab("Empirical Rate") +
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)) 
ggsave("./figures/asym_log_calibration_high.pdf", 
       asym_log_calibration_high,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/asym_log_calibration_high.pdf")

asym_log_high_med <- lapply(asym_log_high_list[1,], function(x) median(x$dep)) |> 
  unlist() |> 
  as_tibble()
asym_log_high_med_boxplot <- asym_log_high_med |> ggplot(aes(y = value, x = 0)) + 
  geom_boxplot(fill = "grey", width = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)) +
  geom_hline(yintercept = 0.1, col = "blue", linetype = 2) +
  ylim(0, 1) +
  xlim(-0.5, 0.5)
ggsave("./figures/asym_log_high_med_boxplot.pdf", 
       asym_log_high_med_boxplot,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/asym_log_high_med_boxplot.pdf")

asym_log_mid_med <- lapply(asym_log_mid_list[1,], function(x) median(x$dep)) |> 
  unlist() |> 
  as_tibble()
asym_log_mid_med_boxplot <- asym_log_mid_med |> ggplot(aes(y = value, x = 0)) + 
  geom_boxplot(fill = "grey", width = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)) +
  geom_hline(yintercept = 0.5, col = "blue", linetype = 2) +
  ylim(0, 1) +
  xlim(-0.5, 0.5)
ggsave("./figures/asym_log_mid_med_boxplot.pdf", 
       asym_log_mid_med_boxplot,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/asym_log_mid_med_boxplot.pdf")

asym_log_low_med <- lapply(asym_log_low_list[1,], function(x) median(x$dep)) |> 
  unlist() |> 
  as_tibble()
asym_log_low_med_boxplot <- asym_log_low_med |> ggplot(aes(y = value, x = 0)) + 
  geom_boxplot(fill = "grey", width = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)) +
  geom_hline(yintercept = 0.9, col = "blue", linetype = 2) +
  ylim(0, 1) +
  xlim(-0.5, 0.5)
ggsave("./figures/asym_log_low_med_boxplot.pdf", 
       asym_log_low_med_boxplot,
       bg = "transparent",
       dpi = 320)
knitr::plot_crop("./figures/asym_log_low_med_boxplot.pdf")

