library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(tidyverse)
library(loo)
library(patchwork)
library(qs)
library(RcppSimdJson)
library(gaugeDependence)
library(ggnewscale)
# extract_posterior_est <- function(dep_type, dep_level, gauge, threshold, likelihood, data_num) {
#   filepath <- sprintf("fits_and_weights/post_params_joint/%s_%s_%s_%s_radial.qs",
#                       gauge, dep_type, dep_level, likelihood, )
#   csvfiles <- paste0(filepath,
#                      list.files(path = filepath, 
#                                 pattern = paste0(dep_level, "_", data_num, "_", likelihood, "_", threshold, "_\\d{1}.csv")))
#   if (gauge != "dirichlet") {
#     temp <- read_cmdstan_csv(csvfiles, variables = "dep")$post_warmup_draws |> as_draws_df()
#     return(median(temp$dep))
#   } else {
#     temp <- read_cmdstan_csv(csvfiles, variables = c("theta1", "theta2"))$post_warmup_draws |> as_draws_df()
#     return(apply(temp[1:2], 2, median))
#   }
# }

## One dataset's gauge fits.  --------
data_num <- 36
mock_data <- fload(sprintf("data/gauss/mid_%s.json", data_num))
r <- mock_data$R
w <- mock_data$W
r0_w <- mock_data$r0_w
rw_df <- cbind(r, w, r0_w) |> as_tibble() |> mutate(high = as.factor(r > r0_w))

gauss_dep_cens <- qread("fits_and_weights/post_params_joint/gauss_gauss_mid_cens_radial.qs")$dep[data_num]
logistic_dep_cens <- qread("fits_and_weights/post_params_joint/logistic_gauss_mid_cens_radial.qs")$dep[data_num]
inv_log_dep_cens <- qread("fits_and_weights/post_params_joint/inv_log_gauss_mid_cens_radial.qs")$dep[data_num]
asym_log_dep_cens <- qread("fits_and_weights/post_params_joint/asym_log_gauss_mid_cens_radial.qs")$dep[data_num]
dirichlet_dep_cens <- as.numeric(qread("fits_and_weights/post_params_joint/dirichlet_gauss_mid_cens_radial.qs")[data_num, c("theta1", "theta2")])
rectangular_dep_cens <- qread("fits_and_weights/post_params_joint/rectangular_gauss_mid_cens_radial.qs")$dep[data_num]

w_sim <- seq(0,1, length.out = nrow(rw_df))
gw_gauss_cens <- gauss_gauge(w_sim, 1-w_sim, gauss_dep_cens)
gw_logistic_cens <- logistic_gauge(w_sim, 1-w_sim, logistic_dep_cens)
gw_inv_log_cens <- inv_log_gauge(w_sim, 1-w_sim, inv_log_dep_cens)
gw_asym_log_cens <- asym_log_gauge(w_sim, 1-w_sim, asym_log_dep_cens)
gw_dirichlet_cens <- dirichlet_gauge(w_sim, 1-w_sim, dirichlet_dep_cens)
gw_rectangular_cens <- rectangular_gauge(w_sim, 1-w_sim, rectangular_dep_cens)

all_fits <- rw_df |> as_tibble() |> cbind(w_sim, gw_gauss_cens, gw_logistic_cens, gw_inv_log_cens, gw_asym_log_cens, gw_dirichlet_cens, gw_rectangular_cens) |> 
  pivot_longer(cols = 6:11, names_to = "gauge_fit", values_to = "values") |>
  mutate(gauge_fit = case_when(gauge_fit == 'gw_logistic_cens' ~ 'Logistic',
                               gauge_fit == 'gw_gauss_cens' ~ 'Gaussian',
                               gauge_fit == 'gw_inv_log_cens' ~ 'Inv. logistic',
                               gauge_fit == 'gw_asym_log_cens' ~ 'Asym. logistic',
                               gauge_fit == 'gw_dirichlet_cens' ~ 'Dirichlet',
                               gauge_fit == 'gw_rectangular_cens' ~ 'Rectangular',
                               .default = gauge_fit))
cols <- c("lightblue", "blue")

rw_fits_cens <- all_fits |> ggplot(aes(x = w, y = r/log(10000), color = high)) + 
  
  # points above and below the threshold
  geom_point(alpha=0.8, aes(color = high)) +
  scale_color_manual(values = cols, guide = "none") +
  
  new_scale_color() + # need this to use a different scale for the gauge function fits
  
  # fits of all 6 gauge functions
  geom_path(aes(x = w_sim, y = 1/values, group = gauge_fit, color = gauge_fit), linewidth = 1) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  
  theme_classic() +
  scale_x_continuous(limits = c(0, 1.01), expand = expansion(mult = c(0,0))) + 
  scale_y_continuous(limits = c(0, 2.01), expand = expansion(mult = c(0,0))) + 
  theme(legend.position = c(0.9, 0.9),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        legend.background = element_rect(fill='transparent', color='transparent'),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  xlab(expression("W"["1"])) + ylab("R/log(n)") + labs(color = "")


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
