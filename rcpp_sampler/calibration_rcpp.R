library(tidyverse)
library(patchwork)
Rcpp::sourceCpp("rcpp_sampler/gauge_mcmc.cpp")

quantile_df <- function(x, probs = c(0.25, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = c('lower', 'upper')
  ) |> pivot_wider(names_from = quant, values_from = val)
}

create_fit_df <- function(gauge, dep_level, lhood_type, data_num) {
  burnin <- ifelse(lhood_type == "trunc", 1000, 500)
  filepath <- paste0("rcpp_sampler/csv_fits/", gauge, "/",
                     list.files(path = paste0("rcpp_sampler/csv_fits/", gauge, "/"), 
                                pattern = paste0(lhood_type, "_", dep_level, "_", data_num, ".csv")))
  return(read_csv(filepath, show_col_types = FALSE) |> as_tibble() |> filter(iter > burnin))
}

create_tib_ci <- function(fit_df, ci_level, truth = 0.5) {
  alpha <- (1-ci_level)/2
  probs_vec <- c(alpha, 1 - alpha)
  temp <- quantile_df(fit_df$theta, probs = probs_vec)
  return(temp |> mutate(level = ci_level, 
                        truth = truth, 
                        coverage = ifelse(lower <= truth & upper >= truth, 1, 0)))
}

coverage_list <- function(gauge, dep_level, lhood_type, data_num, truth) {
  fit_df <- create_fit_df(gauge, dep_level, lhood_type, data_num)
  levels_vec <- seq(0.05, 0.95, by = 0.05)
  coverage_list <- sapply(levels_vec, function(x) create_tib_ci(fit_df, x, truth), simplify = FALSE) |> 
    bind_rows() |>
    mutate(dataset = data_num)
  return(list(fit_df = fit_df, coverage_list = coverage_list))
}

extract_coverage <- function(gauge, dep_level, lhood_type, true_val = 0.5) {
  temp_list <- sapply(1:100, function(y) coverage_list(gauge, dep_level, lhood_type, y, true_val))
  cov_plot_df <- temp_list[2,] |>
    bind_rows() |>
    group_by(level, truth) |>
    summarize(p_hat = mean(coverage),
              sd = sd(coverage)) |>
    ungroup() |>
    mutate(se = sd/10,
           lb = p_hat - qnorm(0.975) * se,
           ub = p_hat + qnorm(0.975) * se,
           truth = as.factor(truth))
  med_mean_df <- lapply(temp_list[1,], function(y) c("median" = median(y$theta), "mean" = mean(y$theta))) |> 
    bind_rows() |> 
    mutate(truth = true_val)
  return(list(cov_plot_df = cov_plot_df,
              med_mean_df = med_mean_df))
}

## Extract coverage from fits -----------
gauss_trunc_mid <- extract_coverage("gauss", "mid", "trunc")
gauss_cens_mid <- extract_coverage("gauss", "mid", "cens")

logistic_trunc_mid <- extract_coverage("logistic", "mid", "trunc")
logistic_cens_mid <- extract_coverage("logistic", "mid", "cens")


## Create boxplots and coverage plots ----------
plot_coverage <- function(cov_tibble) {
  cov_tibble |> ggplot(aes(x = level, y = p_hat)) + 
    geom_point() +
    geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
    geom_abline() +
    theme_classic() +
    xlab("Nominal Rate") +
    ylab("Empirical Rate") +
    scale_x_continuous(limit = c(-0.01, 1.01), expand = c(0,0)) +
    scale_y_continuous(limit = c(-0.01, 1.01), expand = c(0,0)) +
    theme(panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.text = element_text(size = rel(1.2)),
          axis.title = element_text(size = rel(1.2))) 
}

plot_boxplot <- function(cov_tibble, true_val) {
  cov_tibble |> filter(truth == true_val) |>
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
    geom_hline(yintercept = true_val, col = "blue", linetype = 2, linewidth = 1.2) +
    ylim(0, 1) +
    xlim(-0.5, 0.5)
}

column_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = 'Gaussian',
                                                          face = "bold",
                                                          size = 15))
column_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = 'Logistic',
                                                          face = "bold",
                                                          size = 15))
trunc_plots <- ((column_label_1 | column_label_2) + plot_layout(widths = c(.95, .95))) /
  ((plot_coverage(gauss_trunc_mid[[1]]) | plot_coverage(logistic_trunc_mid[[1]])) + plot_layout(widths = c(.95, .95))) /
  ((plot_boxplot(gauss_trunc_mid[[2]], 0.5) | plot_boxplot(logistic_trunc_mid[[2]], 0.5)) + plot_layout(widths = c(.95, .95))) /
  plot_layout(heights = c(0.1, 0.95, 0.95)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("~/Desktop/csu/classes/2024/600_koslovsky/project/trunc_calibration.png",
       plot = trunc_plots,
       dpi = 320,
       bg = "transparent",
       width = 8.1, height = 8)
knitr::plot_crop("~/Desktop/csu/classes/2024/600_koslovsky/project/trunc_calibration.png")

cens_plots <- ((column_label_1 | column_label_2) + plot_layout(widths = c(.95, .95))) /
  ((plot_coverage(gauss_cens_mid[[1]]) | plot_coverage(logistic_cens_mid[[1]])) + plot_layout(widths = c(.95, .95))) /
  ((plot_boxplot(gauss_cens_mid[[2]], 0.5) | plot_boxplot(logistic_cens_mid[[2]], 0.5)) + plot_layout(widths = c(.95, .95))) /
  plot_layout(heights = c(0.1, 0.95, 0.95)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("~/Desktop/csu/classes/2024/600_koslovsky/project/cens_calibration.png",
       plot = cens_plots,
       dpi = 320,
       bg = "transparent",
       width = 8.1, height = 8)
knitr::plot_crop("~/Desktop/csu/classes/2024/600_koslovsky/project/cens_calibration.png")

## Real gauge function with 100 gauge fits to median of dependnece parameters---------
w <- seq(0, 1, length.out = 500)
gw_gauss_truth <- gauss_gauge(w, 0.5)
w_gw_gauss_truth <- cbind(w, gw_gauss_truth)
colnames(w_gw_gauss_truth) <- c("w", "gw")

gw_gauss_100fits_trunc <- sapply(gauss_trunc_mid[[2]]$median, function(x) gauss_gauge(w, x)) |> 
  as_tibble() |> 
  mutate(w1 = w, w2 = 1 - w) |> 
  pivot_longer(cols = !c(w1, w2), names_to = "dataset", values_to = "gw") |>
  mutate(dataset = as.factor(gsub("V", "", dataset)))
gw_gauss_100fits_trunc_plot <- gw_gauss_100fits_trunc |>
  ggplot(aes(x = w1/gw, y = w2/gw, color = dataset)) + geom_path(alpha = 0.5) +
  geom_path(data = w_gw_gauss_truth, 
            mapping = aes(x = w/gw, y = (1-w)/gw), inherit.aes = FALSE, linewidth = 0.75) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = "none") + 
  xlab(expression("W"["1"])) + 
  ylab(expression("W"["2"]))

gw_gauss_100fits_cens <- sapply(gauss_cens_mid[[2]]$median, function(x) gauss_gauge(w, x)) |> 
  as_tibble() |> 
  mutate(w1 = w, w2 = 1 - w) |> 
  pivot_longer(cols = !c(w1, w2), names_to = "dataset", values_to = "gw") |>
  mutate(dataset = as.factor(gsub("V", "", dataset)))
gw_gauss_100fits_cens_plot <- gw_gauss_100fits_cens |>
  ggplot(aes(x = w1/gw, y = w2/gw, color = dataset)) + geom_path(alpha = 0.5) +
  geom_path(data = w_gw_gauss_truth, 
            mapping = aes(x = w/gw, y = (1-w)/gw), inherit.aes = FALSE, linewidth = 0.75) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = "none") + 
  xlab(expression("W"["1"])) + 
  ylab(expression("W"["2"]))

gw_logistic_truth <- logistic_gauge(w, 0.5)
w_gw_logistic_truth <- cbind(w, gw_logistic_truth)
colnames(w_gw_logistic_truth) <- c("w", "gw")

gw_logistic_100fits_trunc <- sapply(logistic_trunc_mid[[2]]$median, function(x) logistic_gauge(w, x)) |> 
  as_tibble() |> 
  mutate(w1 = w, w2 = 1 - w) |> 
  pivot_longer(cols = !c(w1, w2), names_to = "dataset", values_to = "gw") |>
  mutate(dataset = as.factor(gsub("V", "", dataset)))
gw_logistic_100fits_trunc_plot <- gw_logistic_100fits_trunc |>
  ggplot(aes(x = w1/gw, y = w2/gw, color = dataset)) + geom_path(alpha = 0.5) +
  geom_path(data = w_gw_logistic_truth, 
            mapping = aes(x = w/gw, y = (1-w)/gw), inherit.aes = FALSE, linewidth = 0.75) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = "none") + 
  xlab(expression("W"["1"])) + 
  ylab(expression("W"["2"]))


gw_logistic_100fits_cens <- sapply(logistic_cens_mid[[2]]$median, function(x) logistic_gauge(w, x)) |> 
  as_tibble() |> 
  mutate(w1 = w, w2 = 1 - w) |> 
  pivot_longer(cols = !c(w1, w2), names_to = "dataset", values_to = "gw") |>
  mutate(dataset = as.factor(gsub("V", "", dataset)))
gw_logistic_100fits_cens_plot <- gw_logistic_100fits_cens |>
  ggplot(aes(x = w1/gw, y = w2/gw, color = dataset)) + geom_path(alpha = 0.5) +
  geom_path(data = w_gw_logistic_truth, 
            mapping = aes(x = w/gw, y = (1-w)/gw), inherit.aes = FALSE, linewidth = 0.75) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = "none") + 
  xlab(expression("W"["1"])) + 
  ylab(expression("W"["2"]))
  

gauge_plots <- ((column_label_1 | column_label_2) + plot_layout(widths = c(.95, .95))) /
  ((gw_gauss_100fits_trunc_plot | gw_logistic_100fits_trunc_plot) + plot_layout(widths = c(.95, .95))) /
  ((gw_gauss_100fits_cens_plot | gw_logistic_100fits_cens_plot) + plot_layout(widths = c(.95, .95))) +
  plot_layout(heights = c(0.1, 0.95, 0.95)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("~/Desktop/csu/classes/2024/600_koslovsky/project/gauge_fits.png",
       plot = gauge_plots,
       dpi = 320,
       bg = "transparent",
       width = 8.1, height = 8)
knitr::plot_crop("~/Desktop/csu/classes/2024/600_koslovsky/project/trunc_calibration.png")

gauss_trunc_fit <- create_fit_df("gauss", "mid", "trunc", sample(1:100, 1))
gauss_trunc_fit |> ggplot(aes(x = iter, y = alpha, color = as.factor(chain))) + geom_line(alpha = 0.5)

create_mcmc_fit <- function(gauge, dep_level, lhood_type, data_num) {
  filepath <- paste0("rcpp_sampler/csv_fits/", gauge, "/",
                     list.files(path = paste0("rcpp_sampler/csv_fits/", gauge, "/"), 
                                pattern = paste0(lhood_type, "_", dep_level, "_", data_num, ".csv")))
  return(read_csv(filepath, show_col_types = FALSE) |> as_tibble())
}

gauss_trunc_fit <- create_mcmc_fit("gauss", "mid", "trunc", 5)
alpha_gauss_trunc_fit <- gauss_trunc_fit |>
  ggplot(aes(x = iter, y = alpha, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")

theta_gauss_trunc_fit <- gauss_trunc_fit |>
  ggplot(aes(x = iter, y = theta, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")

gauss_cens_fit <- create_mcmc_fit("gauss", "mid", "cens", 5)
alpha_gauss_cens_fit <-  gauss_cens_fit |> 
  ggplot(aes(x = iter, y = alpha, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")

theta_gauss_cens_fit <-  gauss_cens_fit |>  
  ggplot(aes(x = iter, y = theta, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")

logistic_trunc_fit <- create_mcmc_fit("logistic", "mid", "trunc", 5)
alpha_logistic_trunc_fit <- logistic_trunc_fit |>
  ggplot(aes(x = iter, y = alpha, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")

theta_logistic_trunc_fit <- logistic_trunc_fit |>
  ggplot(aes(x = iter, y = theta, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")

logistic_cens_fit <- create_mcmc_fit("logistic", "mid", "cens", 5)
alpha_logistic_cens_fit <-  logistic_cens_fit |> 
  ggplot(aes(x = iter, y = alpha, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")

theta_logistic_cens_fit <-  logistic_cens_fit |>  
  ggplot(aes(x = iter, y = theta, color = as.factor(chain))) + 
  geom_line() + 
  theme_classic() +
  theme(legend.position = "none") + 
  xlab("Iteration") + 
  ylab("Parameter Estimate")


library(patchwork)
column_label_alpha <- wrap_elements(panel = ggpubr::text_grob(label = expression(alpha),
                                                                 size = 15))
column_label_theta <- wrap_elements(panel = ggpubr::text_grob(label = expression(theta),
                                                                 size = 15))
row_label_gauss <- wrap_elements(panel = ggpubr::text_grob(label = "Gaussian",
                                                           size = 15))
row_label_logistic <- wrap_elements(panel = ggpubr::text_grob(label = "Logistic",
                                                           size = 15))
trunc_traceplots <- ((plot_spacer() | column_label_alpha | column_label_theta) + plot_layout(widths = c(0.25, 0.95, 0.95))) / 
  ((row_label_gauss | alpha_gauss_trunc_fit | theta_gauss_trunc_fit) + plot_layout(widths = c(0.25, 0.95, 0.95))) /
  ((row_label_logistic | alpha_logistic_trunc_fit | theta_logistic_trunc_fit) + plot_layout(widths = c(0.25, 0.95, 0.95))) + 
  plot_layout(heights = c(0.15, 1, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/classes/2024/600_koslovsky/project/trunc_traceplots.png",
       plot = trunc_traceplots,
       dpi = 320,
       bg = "transparent",
       width = 10, height = 9)
knitr::plot_crop("~/Desktop/csu/classes/2024/600_koslovsky/project/trunc_traceplots.png")

cens_traceplots <- ((plot_spacer() | column_label_alpha | column_label_theta) + plot_layout(widths = c(0.25, 0.95, 0.95))) / 
  ((row_label_gauss | alpha_gauss_cens_fit | theta_gauss_cens_fit) + plot_layout(widths = c(0.25, 0.95, 0.95))) /
  ((row_label_logistic | alpha_logistic_cens_fit | theta_logistic_cens_fit) + plot_layout(widths = c(0.25, 0.95, 0.95))) + 
  plot_layout(heights = c(0.15, 1, 1)) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("~/Desktop/csu/classes/2024/600_koslovsky/project/cens_traceplots.png",
       plot = cens_traceplots,
       dpi = 320,
       bg = "transparent",
       width = 10, height = 9)
knitr::plot_crop("~/Desktop/csu/classes/2024/600_koslovsky/project/cens_traceplots.png")

library(coda)
gauss_trunc_chains <- split(gauss_trunc_fit, gauss_trunc_fit$chain)
gauss_trunc_mcmc_objects <- lapply(gauss_trunc_chains, function(x) {
  # Remove the 'chain' and 'iteration' column for the mcmc object, assume iterations are in order
  mcmc(as.matrix(x[, c("alpha", "theta")]), start = 1, thin = 1)
})
gauss_trunc_mcmc_list <- mcmc.list(gauss_trunc_mcmc_objects)
round(autocorr.diag(gauss_trunc_mcmc_list), 3)

gauss_cens_chains <- split(gauss_cens_fit, gauss_cens_fit$chain)
gauss_cens_mcmc_objects <- lapply(gauss_cens_chains, function(x) {
  # Remove the 'chain' and 'iteration' column for the mcmc object, assume iterations are in order
  mcmc(as.matrix(x[, c("alpha", "theta")]), start = 1, thin = 1)
})
gauss_cens_mcmc_list <- mcmc.list(gauss_cens_mcmc_objects)
round(autocorr.diag(gauss_cens_mcmc_list), 3)

logistic_trunc_chains <- split(logistic_trunc_fit, logistic_trunc_fit$chain)
logistic_trunc_mcmc_objects <- lapply(logistic_trunc_chains, function(x) {
  # Remove the 'chain' and 'iteration' column for the mcmc object, assume iterations are in order
  mcmc(as.matrix(x[, c("alpha", "theta")]), start = 1, thin = 1)
})
logistic_trunc_mcmc_list <- mcmc.list(logistic_trunc_mcmc_objects)
round(autocorr.diag(logistic_trunc_mcmc_list), 3)

logistic_cens_chains <- split(logistic_cens_fit, logistic_cens_fit$chain)
logistic_cens_mcmc_objects <- lapply(logistic_cens_chains, function(x) {
  # Remove the 'chain' and 'iteration' column for the mcmc object, assume iterations are in order
  mcmc(as.matrix(x[, c("alpha", "theta")]), start = 1, thin = 1)
})
logistic_cens_mcmc_list <- mcmc.list(logistic_cens_mcmc_objects)
round(autocorr.diag(logistic_cens_mcmc_list), 3)
