library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)
library(patchwork)
library(progressr)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

## calibration boxplots -------
extract_param_est <- function(gauge, dep_level, lhood_type, thresh_type, data_num) {
  csvfiles <- list.files(path = paste0("stan/radial_angular/csv_fits/calibrate/", gauge, "/"), 
                         pattern = paste0(dep_level, "_", data_num, "_", lhood_type, "_", thresh_type, "_\\d{1}.csv"),
                         full.names = TRUE)
  fit <- as_cmdstan_fit(csvfiles)
  if (gauge != "dirichlet") {
    param <- fit |> as_draws_df() |>
      select(-c(".iteration", ".chain")) |>
      rename(draw = ".draw") |> 
      select(dep)
    return(tibble(param_est = median(param$dep),
                  dataset_num = data_num))
  } else {
    params <- fit |> as_draws_df() |>
      select(-c(".iteration", ".chain")) |>
      rename(draw = ".draw") |> 
      select(theta1, theta2)
    return(apply(params, 2, function(col) median(col)) |> 
             as_tibble() |> 
             mutate(dataset_num = data_num))
  }
}

# combos <- expand_grid(thres = c("marg", "ctau"),
#                       lhood = c("trunc", "cens"))
# test <- apply(combos, 1, function(row) {
#   sapply(1:3, function(x) extract_param_est("gauss", "high", row["lhood"], row["thres"], x)) |>
#     t() |>
#     as_tibble() |>
#     mutate(param_est = unlist(param_est),
#            dataset_num = unlist(dataset_num),
#            truth = as.numeric(test_truth[test_level]),
#            model_type = paste0(row["lhood"], ", ",  row["thres"]))
# })

create_tib_ests <- function(gauge, dep_level) {
  if(gauge == "gauss") {
    true_vals <- list(high = 0.9, mid = 0.5, low = 0.1)
  } else {
    true_vals <- list(high = 0.1, mid = 0.5, low = 0.9)
  }
  lhood_thres_combo <- expand_grid(thres = c("marg", "ctau"),
                                   lhood = c("trunc", "cens"))
  one_level_ests <- apply(lhood_thres_combo, 1, function(row) {
    sapply(1:100, function(x) extract_param_est(gauge, dep_level, 
                                                row["lhood"], row["thres"], x)) |>
      t() |>
      as_tibble() |> 
      mutate(param_est = unlist(param_est),
             dataset_num = unlist(dataset_num),
             truth = as.numeric(true_vals[dep_level]),
             model_type = paste0(row["lhood"], ", ", row["thres"]))
  }) |> 
    bind_rows() |> 
    mutate(model_type = as.factor(model_type))
  return(one_level_ests)
}

make_boxplot <- function(gauge, dep_level) {
  temp <- create_tib_ests(gauge, dep_level)
  temp_plot <- temp |> ggplot(aes(y = param_est, x = model_type, fill = model_type)) + 
    geom_boxplot() + 
    ylim(0, 1) + 
    xlab("Model") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = rel(1.0)),
          axis.text.x = element_text(size = rel(0.9)),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.text = element_text(size = rel(1.2)),
          axis.title = element_text(size = rel(1.2))) +
    geom_hline(aes(yintercept = truth), col = "darkgrey", linetype = 2, linewidth = 1.2)
  file_name <- paste0(gauge, "_", dep_level,"_boxplot.RDS")
  saveRDS(temp_plot, file = paste0("bma_update_deck/", file_name))
}

dep_level_combos <- expand_grid(dep_type = c("gauss", "logistic"),
                                levels = c("high", "mid", "low"))
with_progress({
  # Create a progress handler
  p <- progressor(steps = nrow(dep_level_combos))
  
  # Apply the function using apply and update the progress bar
  apply(dep_level_combos, 1, function(row) {
    p()  # Update the progress bar
    make_boxplot(gauge = row["dep_type"],
                 dep_level = row["levels"])
  })
})

# now put the dependence levels together in one plot for each dependence type
row_label_1 <- wrap_elements(panel = ggpubr::text_grob(label = "Gaussian",
                                                       face = "bold",
                                                       size = 10))
row_label_2 <- wrap_elements(panel = ggpubr::text_grob(label = "Logistic",
                                                       face = "bold",
                                                       size = 10))
blank_space <- wrap_elements(panel = ggpubr::text_grob(label = '',
                                                       size = 5))

gauss_low_boxplot <- readRDS("~/Desktop/research/gaugeDependence/bma_update_deck/gauss_low_boxplot.RDS")
gauss_mid_boxplot <- readRDS("~/Desktop/research/gaugeDependence/bma_update_deck/gauss_mid_boxplot.RDS")
gauss_high_boxplot <- readRDS("~/Desktop/research/gaugeDependence/bma_update_deck/gauss_high_boxplot.RDS")

logistic_low_boxplot <- readRDS("~/Desktop/research/gaugeDependence/bma_update_deck/logistic_low_boxplot.RDS")
logistic_mid_boxplot <- readRDS("~/Desktop/research/gaugeDependence/bma_update_deck/logistic_mid_boxplot.RDS")
logistic_high_boxplot <- readRDS("~/Desktop/research/gaugeDependence/bma_update_deck/logistic_high_boxplot.RDS")

all_plots <- (((row_label_1 | (gauss_low_boxplot + coord_cartesian(ylim = c(0, 0.5))) | 
                (gauss_mid_boxplot + coord_cartesian(ylim = c(0.25, 0.75))) | 
                (gauss_high_boxplot + coord_cartesian(ylim = c(0.5, 1.0)))) + plot_layout(widths = c(0.15, .95, .95, .95))) /
  ((row_label_2 | (logistic_high_boxplot + coord_cartesian(ylim = c(0, 0.5))) | 
     (logistic_mid_boxplot + coord_cartesian(ylim = c(0.25, 0.75))) | 
     (logistic_low_boxplot + coord_cartesian(ylim = c(0.5, 1.0)))) + plot_layout(widths = c(0.15, .95, .95, .95)))) +
  plot_layout(heights = c(0.85, 0.85)) &
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))

ggsave("bma_update_deck/calibration_boxplot.pdf",
       plot = all_plots,
       dpi = 320,
       bg = "transparent",
       width = 15, height = 8)


## coverage probabilities --------------------
# # following function is for creating coverage probability plots
# quantile_df <- function(x, probs = c(0.25, 0.75)) {
#   tibble(
#     val = quantile(x, probs, na.rm = TRUE),
#     quant = c('lower', 'upper')
#   ) |> pivot_wider(names_from = quant, values_from = val)
# }

# csvfiles <- paste0("stan/csv_fits/calibrate/", "gauss", "/",
#                    list.files(path = paste0("stan/csv_fits/calibrate/", "gauss", "/"), 
#                               pattern = paste0("high", "_", 50, "_", "cens", "_", "marg", "_\\d{1}.csv")))
# fit <- as_cmdstan_fit(csvfiles)
# 
# test <- fit$draws(variables = c("alpha", "dep")) |> as_draws_df() |> select(-c(".iteration", ".chain"))
# plot(test$alpha, test$dep)
# coda::HPDinterval(coda::as.mcmc(test$alpha))
# quantile(test$alpha, c(0.025, 0.975))


# create_fit_df <- function(gauge, dep_level, lhood_type, thresh_type, data_num) {
#   csvfiles <- list.files(path = paste0("stan/csv_fits/calibrate/", gauge, "/"), 
#                          pattern = paste0(dep_level, "_", data_num, "_", lhood_type, "_", thresh_type, "_\\d{1}.csv"),
#                          full.names = TRUE)
#   fit <- as_cmdstan_fit(csvfiles)
#   if (gauge != "dirichlet") {
#     return(fit |> as_draws_df() |>
#              select(-c(".iteration", ".chain")) |>
#              rename(draw = ".draw") |> 
#              select(dep))
#   } else {
#     return(fit |> as_draws_df() |>
#              select(-c(".iteration", ".chain")) |>
#              rename(draw = ".draw") |> 
#              select(theta1, theta2))
#   }
# }

# create_tib_ci <- function(fit_df, ci_level, truth) {
#   alpha <- (1-ci_level)/2
#   probs_vec <- c(alpha, 1 - alpha)
#   if (ncol(fit_df) == 1) {
#     temp <- quantile_df(fit_df$dep, probs = probs_vec)
#     return(temp |> mutate(level = ci_level, 
#                           truth = truth, 
#                           coverage = ifelse(lower <= truth & upper >= truth, 1, 0)))
#   } else {
#     temp <- apply(fit_df, 2, quantile_df, probs = probs_vec) |> 
#       bind_rows() |> 
#       mutate(param = c("theta1", "theta2"), .before = "lower")
#     return(temp |> mutate(level = rep(ci_level, 2), 
#                           truth = c(truth, 2), 
#                           coverage = ifelse(lower <= truth & upper >= truth, 1, 0)))
#   }
# }
# 
# # for dirichlet
# dep_levels <- list(c(3, "high"), c(1, "mid"), c(0.5, "low"))

# coverage_list <- function(gauge, dep_level, lhood_type, thresh_type, data_num, truth) {
#   fit_df <- create_fit_df(gauge, dep_level, lhood_type, thresh_type, data_num)
#   levels_vec <- seq(0.05, 0.95, by = 0.05)
#   cov_list <- sapply(levels_vec, function(x) create_tib_ci(fit_df, x, truth), simplify = FALSE) |> 
#     bind_rows() |>
#     mutate(dataset = data_num)
#   return(list(fit_df = fit_df, cov_list = cov_list))
# }

# extract_coverage <- function(gauge, lhood_type, thresh_type) {
#   temp_list <- vector("list", 3)
#   names(temp_list) <- c("high", "mid", "low")
#   if (gauge == "gauss") {
#     true_vals <- c(0.9, 0.5, 0.1)
#     for(i in seq_along(true_vals)) {
#       temp_list[[i]] <- sapply(1:100, function(y) coverage_list("gauss", names(temp_list)[i], lhood_type, thresh_type, y, true_vals[i]))
#     }
#   } else if (gauge == "dirichlet") {
#     true_vals <- c(3, 1, 0.5)  
#     for(i in seq_along(true_vals)) {
#       temp_list[[i]] <- sapply(1:100, function(y) coverage_list("dirichlet", names(temp_list)[i], lhood_type, thresh_type, y, true_vals[i]))
#     }
#   } else {
#     true_vals <- c(0.1, 0.5, 0.9)
#     for(i in seq_along(true_vals)) {
#       temp_list[[i]] <- sapply(1:100, function(y) coverage_list(gauge, names(temp_list)[i], lhood_type, thresh_type, y, true_vals[i]))
#     }  
#   }
#   cov_plot_df <- lapply(temp_list, function(x) x[2,] |> bind_rows()) |>
#     bind_rows() |>
#     group_by(level, truth) |>
#     summarize(p_hat = mean(coverage),
#               sd = sd(coverage)) |>
#     ungroup() |>
#     mutate(se = sd/10,
#            lb = p_hat - qnorm(0.975) * se,
#            ub = p_hat + qnorm(0.975) * se,
#            truth = as.factor(truth))
#   med_mean_df <- lapply(temp_list, function(x) lapply(x[1,], function(y)  c("median" = median(y$dep), "mean" = mean(y$dep)))) |> bind_rows() |> mutate(truth = rep(true_vals, times = rep(100,3)))
#   return(list(cov_plot_df = cov_plot_df,
#               med_mean_df = med_mean_df))
# }

## NEED TO CREATE ALTERNATIVE FUNCTION TO EXTRACT MEAN AND MEDIAN FOR DIRICHLET GAUGE FUNCTION ####



## Create boxplots and coverage plots ----------
# plot_coverage <- function(cov_tibble, true_dep) {
#   cov_tibble |> filter(truth == true_dep) |> ggplot(aes(x = level, y = p_hat)) + 
#     geom_point() +
#     geom_linerange(aes(ymin = lb, ymax = ub), linetype = 2, col = "red") +
#     geom_abline() +
#     theme_classic() +
#     ylim(-0.01,1.01) +
#     xlim(-0.01,1.01) +
#     xlab("Nominal Rate") +
#     ylab("Empirical Rate") +
#     theme(panel.background = element_rect(fill='transparent'),
#           plot.background = element_rect(fill='transparent', color=NA),
#           axis.text = element_text(size = rel(1.2)),
#           axis.title = element_text(size = rel(1.2))) 
# }

# plot_boxplot <- function(cov_tibble, true_dep) {
#   cov_tibble |> filter(truth == true_dep) |> 
#     ggplot(aes(y = median, x = 0)) + 
#     geom_boxplot(fill = "grey", width = 0.5) +
#     theme_classic() +
#     theme(
#       axis.title.x = element_blank(),
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       axis.title.y = element_blank(),
#       panel.background = element_rect(fill='transparent'),
#       plot.background = element_rect(fill='transparent', color=NA),
#       axis.text = element_text(size = rel(1.2)),
#       axis.title = element_text(size = rel(1.2))) +
#     geom_hline(yintercept = true_dep, col = "blue", linetype = 2, linewidth = 1.2) +
#     ylim(0, 1) +
#     xlim(-0.5, 0.5)
# }

# create_save_plots <- function(coverage_list) {
#   
# }



