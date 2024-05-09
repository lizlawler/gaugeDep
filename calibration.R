library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)
library(patchwork)

quantile_df <- function(x, probs = c(0.25, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = c('lower', 'upper')
  ) |> pivot_wider(names_from = quant, values_from = val)
}

csvfiles <- paste0("stan/csv_fits/calibrate/", "gauss", "/",
                   list.files(path = paste0("stan/csv_fits/calibrate/", "gauss", "/"), 
                              pattern = paste0("high", "_", 50, "_", "cens", "_", "marg", "_\\d{1}.csv")))
fit <- as_cmdstan_fit(csvfiles)

test <- fit$draws(variables = c("alpha", "dep")) |> as_draws_df() |> select(-c(".iteration", ".chain"))
plot(test$alpha, test$dep)
coda::HPDinterval(coda::as.mcmc(test$alpha))
quantile(test$alpha, c(0.025, 0.975))


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

## NEED TO CREATE ALTERNATIVE FUNCTION TO EXTRACT MEAN AND MEDIAN FOR DIRICHLET GAUGE FUNCTION ####

## Extract coverage from fits -----------
gauss_coverage_marg <- extract_coverage("gauss", "trunc", "marg")
logistic_coverage_marg <- extract_coverage("logistic", "trunc", "marg")

gauss_coverage_ctau <- extract_coverage("gauss", "trunc", "ctau")
logistic_coverage_ctau <- extract_coverage("logistic", "trunc", "ctau")

gauss_coverage_cens <- extract_coverage("gauss", "cens", "marg")
logistic_coverage_cens <- extract_coverage("logistic", "cens", "marg")

inv_log_coverage_marg <- extract_coverage("inv_log", "trunc", "marg")
asym_log_coverage_marg <- extract_coverage("asym_log", "trunc", "marg")

inv_log_coverage_ctau <- extract_coverage("inv_log", "trunc", "ctau")
asym_log_coverage_ctau <- extract_coverage("asym_log", "trunc", "ctau")

inv_log_coverage_cens_marg <- extract_coverage("inv_log", "cens", "marg")

dirichlet_coverage_marg <- extract_coverage("dirichlet", "trunc", "marg")

gauss_coverage_cens_take2 <- extract_coverage("gauss", "cens", "marg_take2")

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
    geom_hline(yintercept = true_dep, col = "blue", linetype = 2, linewidth = 1.2) +
    ylim(0, 1) +
    xlim(-0.5, 0.5)
}

# create_save_plots <- function(coverage_list) {
#   
# }

all_gauss_marg_plots <- (plot_coverage(gauss_coverage_marg[[1]], 0.1) | plot_coverage(gauss_coverage_marg[[1]], 0.5) | plot_coverage(gauss_coverage_marg[[1]], 0.9)) / 
  (plot_boxplot(gauss_coverage_marg[[2]], 0.1) | plot_boxplot(gauss_coverage_marg[[2]], 0.5) | plot_boxplot(gauss_coverage_marg[[2]], 0.9))

ggsave("./figures/gauss_marg_calibration.pdf", 
       all_gauss_marg_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

all_gauss_cens_take2_plots <- (plot_coverage(gauss_coverage_cens_take2[[1]], 0.1) | plot_coverage(gauss_coverage_cens_take2[[1]], 0.5) | plot_coverage(gauss_coverage_cens_take2[[1]], 0.9)) / 
  (plot_boxplot(gauss_coverage_cens_take2[[2]], 0.1) | plot_boxplot(gauss_coverage_cens_take2[[2]], 0.5) | plot_boxplot(gauss_coverage_cens_take2[[2]], 0.9))

ggsave("./figures/gauss_cens_take2_calibration.pdf", 
       all_gauss_cens_take2_plots,
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

all_logistic_cens_plots <- (plot_coverage(logistic_coverage_cens[[1]], 0.1) | plot_coverage(logistic_coverage_cens[[1]], 0.5) | plot_coverage(logistic_coverage_cens[[1]], 0.9)) / 
  (plot_boxplot(logistic_coverage_cens[[2]], 0.1) | plot_boxplot(logistic_coverage_cens[[2]], 0.5) | plot_boxplot(logistic_coverage_cens[[2]], 0.9))

ggsave("./figures/logistic_cens_calibration.pdf", 
       all_logistic_cens_plots,
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

all_inv_log_marg_plots <- (plot_coverage(inv_log_coverage_marg[[1]], 0.1) | plot_coverage(inv_log_coverage_marg[[1]], 0.5) | plot_coverage(inv_log_coverage_marg[[1]], 0.9)) / 
  (plot_boxplot(inv_log_coverage_marg[[2]], 0.1) | plot_boxplot(inv_log_coverage_marg[[2]], 0.5) | plot_boxplot(inv_log_coverage_marg[[2]], 0.9))

ggsave("./figures/inv_log_marg_calibration.pdf", 
       all_inv_log_marg_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

all_inv_log_cens_marg_plots <- (plot_coverage(inv_log_coverage_cens_marg[[1]], 0.1) | plot_coverage(inv_log_coverage_cens_marg[[1]], 0.5) | plot_coverage(inv_log_coverage_cens_marg[[1]], 0.9)) / 
  (plot_boxplot(inv_log_coverage_cens_marg[[2]], 0.1) | plot_boxplot(inv_log_coverage_cens_marg[[2]], 0.5) | plot_boxplot(inv_log_coverage_cens_marg[[2]], 0.9))

ggsave("./figures/inv_log_cens_marg_calibration.pdf", 
       all_inv_log_cens_marg_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

all_asym_log_marg_plots <- (plot_coverage(asym_log_coverage_marg[[1]], 0.1) | plot_coverage(asym_log_coverage_marg[[1]], 0.5) | plot_coverage(asym_log_coverage_marg[[1]], 0.9)) / 
  (plot_boxplot(asym_log_coverage_marg[[2]], 0.1) | plot_boxplot(asym_log_coverage_marg[[2]], 0.5) | plot_boxplot(asym_log_coverage_marg[[2]], 0.9))

ggsave("./figures/asym_log_marg_calibration.pdf", 
       all_asym_log_marg_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)


all_inv_log_ctau_plots <- (plot_coverage(inv_log_coverage_ctau[[1]], 0.1) | plot_coverage(inv_log_coverage_ctau[[1]], 0.5) | plot_coverage(inv_log_coverage_ctau[[1]], 0.9)) / 
  (plot_boxplot(inv_log_coverage_ctau[[2]], 0.1) | plot_boxplot(inv_log_coverage_ctau[[2]], 0.5) | plot_boxplot(inv_log_coverage_ctau[[2]], 0.9))

ggsave("./figures/inv_log_ctau_calibration.pdf", 
       all_inv_log_ctau_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

out <- (plot_coverage(asym_log_coverage_ctau[[1]], 0.1) | plot_coverage(asym_log_coverage_ctau[[1]], 0.5) | plot_coverage(asym_log_coverage_ctau[[1]], 0.9)) / 
  (plot_boxplot(asym_log_coverage_ctau[[2]], 0.1) | plot_boxplot(asym_log_coverage_ctau[[2]], 0.5) | plot_boxplot(asym_log_coverage_ctau[[2]], 0.9))

ggsave("./figures/asym_log_ctau_calibration.pdf", 
       all_asym_log_ctau_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)

all_dirichlet_marg_plots <- (plot_coverage(dirichlet_coverage_marg[[1]], 0.5) | plot_coverage(dirichlet_coverage_marg[[1]], 1) | plot_coverage(dirichlet_coverage_marg[[1]], 2) | plot_coverage(dirichlet_coverage_marg[[1]], 3)) / 
  (plot_boxplot(dirichlet_coverage_marg[[2]], 0.5) | plot_boxplot(dirichlet_coverage_marg[[2]], 1) | plot_boxplot(dirichlet_coverage_marg[[2]], 2) | plot_boxplot(dirichlet_coverage_marg[[2]], 3))

ggsave("./figures/dirichlet_marg_calibration.pdf", 
       all_dirichlet_marg_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)


