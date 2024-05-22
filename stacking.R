library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(tidyverse)
library(loo)
library(patchwork)
library(RcppSimdJson)

## Gauge functions
gauss_gauge <- function(x, y, rho = 0.5) {
  top <- x + y - 2 * rho * sqrt(x * y)
  return(top/(1-rho^2))
}

logistic_gauge <- function(x, y, r = 0.5) {
  r_inv <- 1/r
  return(r_inv * pmax(x, y) + (1-r_inv)*pmin(x,y))
}

inv_log_gauge <- function(x, y, r = 0.5) ((x^(1/r) + y^(1/r))^r)

asym_log_gauge <- function(x, y, r = 0.5) {
  r_inv <- 1/r
  return(pmin((x + y), (r_inv * pmax(x, y) + (1-r_inv)*pmin(x,y))))
}

dirichlet_gauge <- function(x, y, theta1, theta2) {
  return((1 + theta1 + theta2) * pmax(x, y) - (theta1 * x + theta2 * y))
}

rectangular_gauge <- function(x, y, dep) {
  return(pmax((x - y) / dep, (y - x) / dep, (x+ y) / (2 - dep)))
}

extract_posterior_est <- function(dep_type, dep_level, gauge, threshold, likelihood, data_num) {
  filepath <- paste0("stan/csv_fits/stacking/", dep_type, "/", gauge, "/")
  csvfiles <- paste0(filepath,
                     list.files(path = filepath, 
                                pattern = paste0(dep_level, "_", data_num, "_", likelihood, "_", threshold, "_\\d{1}.csv")))
  if (gauge != "dirichlet") {
    temp <- read_cmdstan_csv(csvfiles, variables = "dep")$post_warmup_draws |> as_draws_df()
    return(median(temp$dep))
  } else {
    temp <- read_cmdstan_csv(csvfiles, variables = c("theta1", "theta2"))$post_warmup_draws |> as_draws_df()
    return(apply(temp[1:2], 2, median))
  }
}

## One dataset's gauge fits.  --------
mock_data <- fload("data/gauss/high_36.json")
R <- mock_data$R
W <- mock_data$W
rw_df <- cbind(R, W)
xy_df <- cbind(W*R, R - (W*R)) |> as_tibble() |> rename(X = V1, Y = V2)

gauss_dep <- extract_posterior_est("gauss", "high", "gauss", "marg", "cens", 36)
logistic_dep <- extract_posterior_est("gauss", "high", "logistic", "marg", "cens", 36)
inv_log_dep <- extract_posterior_est("gauss", "high", "inv_log", "marg", "cens", 36)
asym_log_dep <- extract_posterior_est("gauss", "high", "asym_log", "marg", "cens", 36)
dirichlet_dep <- extract_posterior_est("gauss", "high", "dirichlet", "marg", "cens", 36)
rectangular_dep <- extract_posterior_est("gauss", "high", "rectangular", "marg", "cens", 36)

w <- seq(0,1, length.out = nrow(rw_df))
gw_gauss <- gauss_gauge(w, 1-w, gauss_dep)
gw_logistic <- logistic_gauge(w, 1-w, logistic_dep)
gw_inv_log <- inv_log_gauge(w, 1-w, inv_log_dep)
gw_asym_log <- asym_log_gauge(w, 1-w, asym_log_dep)
gw_dirichlet <- dirichlet_gauge(w, 1-w, as.numeric(dirichlet_dep[1]), as.numeric(dirichlet_dep[2]))
gw_rectangular <- rectangular_gauge(w, 1-w, rectangular_dep)

all_fits <- rw_df |> as_tibble() |> cbind(xy_df, w, gw_gauss, gw_logistic, gw_inv_log, gw_asym_log, gw_dirichlet, gw_rectangular) |> 
  pivot_longer(cols = 6:11, names_to = "gauge_fit", values_to = "values") |>
  mutate(gauge_fit = case_when(gauge_fit == 'gw_logistic' ~ 'Logistic',
                               gauge_fit == 'gw_gauss' ~ 'Gaussian',
                               gauge_fit == 'gw_inv_log' ~ 'Inv. logistic',
                               gauge_fit == 'gw_asym_log' ~ 'Asym. logistic',
                               gauge_fit == 'gw_dirichlet' ~ 'Dirichlet',
                               gauge_fit == 'gw_rectangular' ~ 'Rectangular',
                               .default = gauge_fit))
rw_fits <- all_fits |> ggplot(aes(x = W, y = R/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w, y = 1/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.01), expand = c(0,0)) + scale_y_continuous(limits = c(0, 2.01), expand = c(0,0)) + 
  theme(legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab("w") + ylab("r/log(n)") +labs(color = "") +
  scale_color_brewer(palette = "Dark2")

xy_fits <- all_fits |> ggplot(aes(x = X/log(10000), y = Y/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w/values, y = (1-w)/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.21), expand = c(0,0)) + scale_y_continuous(limits = c(0, 1.21), expand = c(0,0)) + 
  theme(legend.position = "none",
    # legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab(expression("X"["1"]/log(n))) + ylab(expression("X"["2"]/log(n))) +labs(color = "") + 
  scale_color_brewer(palette = "Dark2")

gauss_gauge_fits_v2 <- rw_fits + guide_area() + xy_fits + plot_layout(guides = 'collect') +
  plot_layout(widths = c(1, 0.4, 1)) & 
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("figures/gauss_high_cens_36_gauge_fits_xy_rw.pdf",
       gauss_gauge_fits_v2,
       dpi = 320,
       bg = 'transparent',
       width = 14, height = 6.5)
knitr::plot_crop("figures/gauss_high_cens_36_gauge_fits_xy_rw.pdf")

gauss_dep <- extract_posterior_est("gauss", "high", "gauss", "marg", "trunc", 36)
logistic_dep <- extract_posterior_est("gauss", "high", "logistic", "marg", "trunc", 36)
inv_log_dep <- extract_posterior_est("gauss", "high", "inv_log", "marg", "trunc", 36)
asym_log_dep <- extract_posterior_est("gauss", "high", "asym_log", "marg", "trunc", 36)
dirichlet_dep <- extract_posterior_est("gauss", "high", "dirichlet", "marg", "trunc", 36)
rectangular_dep <- extract_posterior_est("gauss", "high", "rectangular", "marg", "trunc", 36)

w <- seq(0,1, length.out = nrow(rw_df))
gw_gauss <- gauss_gauge(w, 1-w, gauss_dep)
gw_logistic <- logistic_gauge(w, 1-w, logistic_dep)
gw_inv_log <- inv_log_gauge(w, 1-w, inv_log_dep)
gw_asym_log <- asym_log_gauge(w, 1-w, asym_log_dep)
gw_dirichlet <- dirichlet_gauge(w, 1-w, as.numeric(dirichlet_dep[1]), as.numeric(dirichlet_dep[2]))
gw_rectangular <- rectangular_gauge(w, 1-w, rectangular_dep)

all_fits <- rw_df |> as_tibble() |> cbind(xy_df, w, gw_gauss, gw_logistic, gw_inv_log, gw_asym_log, gw_dirichlet, gw_rectangular) |> 
  pivot_longer(cols = 6:11, names_to = "gauge_fit", values_to = "values") |>
  mutate(gauge_fit = case_when(gauge_fit == 'gw_logistic' ~ 'Logistic',
                               gauge_fit == 'gw_gauss' ~ 'Gaussian',
                               gauge_fit == 'gw_inv_log' ~ 'Inv. logistic',
                               gauge_fit == 'gw_asym_log' ~ 'Asym. logistic',
                               gauge_fit == 'gw_dirichlet' ~ 'Dirichlet',
                               gauge_fit == 'gw_rectangular' ~ 'Rectangular',
                               .default = gauge_fit))
rw_fits <- all_fits |> ggplot(aes(x = W, y = R/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w, y = 1/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.01), expand = c(0,0)) + scale_y_continuous(limits = c(0, 2.01), expand = c(0,0)) + 
  theme(legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab("w") + ylab("r/log(n)") +labs(color = "") +
  scale_color_brewer(palette = "Dark2")

xy_fits <- all_fits |> ggplot(aes(x = X/log(10000), y = Y/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w/values, y = (1-w)/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.21), expand = c(0,0)) + scale_y_continuous(limits = c(0, 1.21), expand = c(0,0)) + 
  theme(legend.position = "none",
        # legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab(expression("X"["1"]/log(n))) + ylab(expression("X"["2"]/log(n))) +labs(color = "") + 
  scale_color_brewer(palette = "Dark2")

gauss_gauge_fits <- rw_fits + guide_area() + xy_fits + plot_layout(guides = 'collect') +
  plot_layout(widths = c(1, 0.4, 1)) & 
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("figures/gauss_high_trunc_marg_36_gauge_fits.pdf",
       gauss_gauge_fits,
       dpi = 320,
       bg = 'transparent',
       width = 14, height = 6.5)
knitr::plot_crop("figures/gauss_high_trunc_marg_36_gauge_fits.pdf")


## fits on logistic dependence structure
mock_data <- fload("data/logistic/high_36.json")
R <- mock_data$R
W <- mock_data$W
rw_df <- cbind(R, W)
xy_df <- cbind(W*R, R - (W*R)) |> as_tibble() |> rename(X = V1, Y = V2)

gauss_dep <- extract_posterior_est("logistic", "high", "gauss", "marg", "cens", 36)
logistic_dep <- extract_posterior_est("logistic", "high", "logistic", "marg", "cens", 36)
inv_log_dep <- extract_posterior_est("logistic", "high", "inv_log", "marg", "cens", 36)
asym_log_dep <- extract_posterior_est("logistic", "high", "asym_log", "marg", "cens", 36)
dirichlet_dep <- extract_posterior_est("logistic", "high", "dirichlet", "marg", "cens", 36)
rectangular_dep <- extract_posterior_est("logistic", "high", "rectangular", "marg", "cens", 36)

w <- seq(0,1, length.out = nrow(rw_df))
gw_gauss <- gauss_gauge(w, 1-w, gauss_dep)
gw_logistic <- logistic_gauge(w, 1-w, logistic_dep)
gw_inv_log <- inv_log_gauge(w, 1-w, inv_log_dep)
gw_asym_log <- asym_log_gauge(w, 1-w, asym_log_dep)
gw_dirichlet <- dirichlet_gauge(w, 1-w, as.numeric(dirichlet_dep[1]), as.numeric(dirichlet_dep[2]))
gw_rectangular <- rectangular_gauge(w, 1-w, rectangular_dep)

all_fits <- rw_df |> as_tibble() |> cbind(xy_df, w, gw_gauss, gw_logistic, gw_inv_log, gw_asym_log, gw_dirichlet, gw_rectangular) |> 
  pivot_longer(cols = 6:11, names_to = "gauge_fit", values_to = "values") |>
  mutate(gauge_fit = case_when(gauge_fit == 'gw_logistic' ~ 'Logistic',
                               gauge_fit == 'gw_gauss' ~ 'Gaussian',
                               gauge_fit == 'gw_inv_log' ~ 'Inv. logistic',
                               gauge_fit == 'gw_asym_log' ~ 'Asym. logistic',
                               gauge_fit == 'gw_dirichlet' ~ 'Dirichlet',
                               gauge_fit == 'gw_rectangular' ~ 'Rectangular',
                               .default = gauge_fit))
rw_fits <- all_fits |> ggplot(aes(x = W, y = R/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w, y = 1/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.01), expand = c(0,0)) + scale_y_continuous(limits = c(0, 2.01), expand = c(0,0)) + 
  theme(legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab("w") + ylab("r/log(n)") +labs(color = "") +
  scale_color_brewer(palette = "Dark2")

xy_fits <- all_fits |> ggplot(aes(x = X/log(10000), y = Y/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w/values, y = (1-w)/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.21), expand = c(0,0)) + scale_y_continuous(limits = c(0, 1.21), expand = c(0,0)) + 
  theme(legend.position = "none",
        # legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab(expression("X"["1"]/log(n))) + ylab(expression("X"["2"]/log(n))) +labs(color = "") + 
  scale_color_brewer(palette = "Dark2")

logistic_gauge_fits <- rw_fits + guide_area() + xy_fits + plot_layout(guides = 'collect') +
  plot_layout(widths = c(1, 0.4, 1)) & 
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("figures/logistic_high_cens_marg_36_gauge_fits.pdf",
       logistic_gauge_fits,
       dpi = 320,
       bg = 'transparent',
       width = 14, height = 6.5)
knitr::plot_crop("figures/logistic_high_cens_marg_36_gauge_fits.pdf")

gauss_dep <- extract_posterior_est("logistic", "high", "gauss", "marg", "trunc", 36)
logistic_dep <- extract_posterior_est("logistic", "high", "logistic", "marg", "trunc", 36)
inv_log_dep <- extract_posterior_est("logistic", "high", "inv_log", "marg", "trunc", 36)
asym_log_dep <- extract_posterior_est("logistic", "high", "asym_log", "marg", "trunc", 36)
dirichlet_dep <- extract_posterior_est("logistic", "high", "dirichlet", "marg", "trunc", 36)
rectangular_dep <- extract_posterior_est("logistic", "high", "rectangular", "marg", "trunc", 36)

w <- seq(0,1, length.out = nrow(rw_df))
gw_gauss <- gauss_gauge(w, 1-w, gauss_dep)
gw_logistic <- logistic_gauge(w, 1-w, logistic_dep)
gw_inv_log <- inv_log_gauge(w, 1-w, inv_log_dep)
gw_asym_log <- asym_log_gauge(w, 1-w, asym_log_dep)
gw_dirichlet <- dirichlet_gauge(w, 1-w, as.numeric(dirichlet_dep[1]), as.numeric(dirichlet_dep[2]))
gw_rectangular <- rectangular_gauge(w, 1-w, rectangular_dep)

all_fits <- rw_df |> as_tibble() |> cbind(xy_df, w, gw_gauss, gw_logistic, gw_inv_log, gw_asym_log, gw_dirichlet, gw_rectangular) |> 
  pivot_longer(cols = 6:11, names_to = "gauge_fit", values_to = "values") |>
  mutate(gauge_fit = case_when(gauge_fit == 'gw_logistic' ~ 'Logistic',
                               gauge_fit == 'gw_gauss' ~ 'Gaussian',
                               gauge_fit == 'gw_inv_log' ~ 'Inv. logistic',
                               gauge_fit == 'gw_asym_log' ~ 'Asym. logistic',
                               gauge_fit == 'gw_dirichlet' ~ 'Dirichlet',
                               gauge_fit == 'gw_rectangular' ~ 'Rectangular',
                               .default = gauge_fit))
rw_fits <- all_fits |> ggplot(aes(x = W, y = R/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w, y = 1/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.01), expand = c(0,0)) + scale_y_continuous(limits = c(0, 2.01), expand = c(0,0)) + 
  theme(legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab("w") + ylab("r/log(n)") +labs(color = "") +
  scale_color_brewer(palette = "Dark2")

xy_fits <- all_fits |> ggplot(aes(x = X/log(10000), y = Y/log(10000))) + geom_point(color="black", size = 1) +
  geom_path(aes(x = w/values, y = (1-w)/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.21), expand = c(0,0)) + scale_y_continuous(limits = c(0, 1.21), expand = c(0,0)) + 
  theme(legend.position = "none",
        # legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.3)),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.3))) +
  xlab(expression("X"["1"]/log(n))) + ylab(expression("X"["2"]/log(n))) +labs(color = "") + 
  scale_color_brewer(palette = "Dark2")

logistic_gauge_fits <- rw_fits + guide_area() + xy_fits + plot_layout(guides = 'collect') +
  plot_layout(widths = c(1, 0.4, 1)) & 
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("figures/logistic_high_trunc_marg_36_gauge_fits.pdf",
       logistic_gauge_fits,
       dpi = 320,
       bg = 'transparent',
       width = 14, height = 6.5)
knitr::plot_crop("figures/logistic_high_trunc_marg_36_gauge_fits.pdf")



# mod1 <- cmdstan_model("stan/bivar_trunc_gamma_asym_log.stan", force_recompile = TRUE)
# fit_dummy <- mod1$sample(data = "data/independent/low_1.json",
#                          iter_warmup = 1,
#                          iter_sampling = 1)
# 
# start_file_path <- paste0("stan/csv_fits/", "stacking", "/", "independent", "/", "asym_log", "/")
# csvfiles <- paste0(start_file_path,
#                    list.files(path = start_file_path, 
#                               pattern = paste0("low", "_", 1, "_\\d{1}.csv")))
# fit_csv <- as_cmdstan_fit(csvfiles)
# fit_dummy$.__enclos_env__$private$draws_ <- fit_csv$.__enclos_env__$private$draws_
# loo1 <- fit_dummy$loo()
# loo2 <- fit_dummy$loo(moment_match = TRUE)

# MCMCvis::MCMCtrace(fit_csv,
#                    params = c('dep', 'alpha'),
#                    ind = TRUE,
#                    open_pdf = TRUE)
# 
# ids <- pareto_k_ids(loo1)
# pareto_k_values(loo1)[ids]
# pareto_k_values(loo2)[ids]


## Create boxplots of stacking weights ----------------
gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
weights_files <- paste0("stacking_weights/", 
                        list.files(path = "stacking_weights/", pattern = ".RDS"))
create_stacking_boxplot <- function(weights_file) {
  plot_filename <- paste0("figures/", 
                     str_remove(basename(weights_file), ".RDS"), ".pdf")
  all_wts <- readRDS(weights_file) |>
    bind_rows() |> 
    mutate(method = rep(gauge_library, 100)) |>
    mutate(stacking = as.numeric(stacking),
           pseudobma_boot = as.numeric(pseudobma_boot),
           pseudobma_noboot = as.numeric(pseudobma_noboot)) |>
    pivot_longer(cols = 1:3, names_to = 'weighting', values_to = 'weights')
  boxplot_wts <- all_wts |>
    mutate(method = case_when(method == "gauss" ~ 'Gauss',
                              method == 'logistic' ~ 'Logistic',
                              method == 'inv_log' ~ 'Inv. logistic',
                              method == 'asym_log' ~ 'Asym. logistic',
                              method == 'rectangular' ~ 'Rectangular',
                              method == 'dirichlet' ~ 'Dirichlet'),
           weighting = case_when(weighting == 'stacking' ~ 'Stacking',
                                 weighting == 'pseudobma_boot' ~ 'Pseudo-BMA+',
                                 weighting == 'pseudobma_noboot' ~ 'Pseudo-BMA'),
           weighting = as.factor(weighting)) |>
    ggplot(aes(x = method, y = weights, fill = weighting)) + geom_boxplot() +
    theme_classic() + 
    # scale_y_continuous(limits = c(0, 1), expand = c(0,0)) + 
    theme(panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'),
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(size = rel(1)),
          axis.text = element_text(size = rel(1)),
          axis.title = element_text(size = rel(1))) +
    xlab("Gauge function") + ylab("Weights") + labs(fill = "")
  ggsave(plot_filename,
         plot = boxplot_wts,
         bg = 'transparent',
         dpi = 320,
         width = 8,
         height = 3.5)
  print(paste0(plot_filename, " has been saved"))
}

sapply(weights_files, create_stacking_boxplot)
