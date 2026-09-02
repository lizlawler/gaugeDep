# =============================================================================
# Assesses calibration of the radial dependence-parameter estimates across the
# 200 simulation datasets. For each gauge/level, produces (1) boxplots of the
# posterior estimates against the true value and (2) empirical-vs-nominal
# coverage plots of the credible intervals, under both likelihoods. Combines
# them side by side per scenario.
#
# Inputs:    fits_and_weights/post_params_joint/{gauge}_{gauge}_{level}_{lhood}_radial.qs
#            samplers/rcpp/radial_mcmc_fits/{gauge}/{gauge}_{lhood}_{level}_{i}.qs
# Outputs:   figures/{gauge}_{level}_coverage.pdf
#            figures/calibration_qs/{gauge}_coverage.qs
# =============================================================================

library(posterior)
library(tidyverse)
library(patchwork)
library(progressr)
library(qs)
library(grafify)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

## calibration boxplots -------
create_tib_ests <- function(gauge, dep_level) {
  if(gauge == "gauss") {
    true_vals <- list(high = 0.9, mid = 0.5, low = 0.1)
  } else {
    true_vals <- list(high = 0.1, mid = 0.5, low = 0.9)
  }
  cens_temp <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_cens_radial.qs", gauge, gauge, dep_level)) |> 
    select(-alpha) |> 
    mutate(lhood = "Censored",
           level = dep_level)
  trunc_temp <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_trunc_radial.qs", gauge, gauge, dep_level)) |> 
    select(-alpha) |> 
    mutate(lhood = "Truncated",
           level = dep_level)
  
  all_ests <- rbind(cens_temp, trunc_temp) |> mutate(truth = as.numeric(true_vals[dep_level]))
  return(all_ests)
}

make_boxplot <- function(gauge, dep_level) {
  temp <- create_tib_ests(gauge, dep_level)
  temp_plot <- temp |> ggplot(aes(y = dep, x = lhood, fill = lhood)) + 
    geom_boxplot() + 
    scale_fill_grafify("r4", ColSeq = FALSE) +
    # scale_x_continuous(expand = expansion(mult = c(0,0.01))) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0,0.01))) +
    xlab("Likelihood") +
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
  # file_name <- paste0(gauge, "_", dep_level,"_boxplot.RDS")
  # saveRDS(temp_plot, file = paste0("bma_update_deck/", file_name))
  obj_name <- sprintf("%s_%s_boxplot", gauge, dep_level)
  assign(obj_name, temp_plot, envir = .GlobalEnv)
  # return(temp_plot)
}

# make_boxplot("gauss", "high")
dep_level_combos <- expand_grid(dep_type = c("gauss", "logistic"),
                                levels = c("high", "mid", "low"))

invisible(apply(dep_level_combos, 1, function(row) {
  make_boxplot(gauge = row["dep_type"],
               dep_level = row["levels"])
}))

## coverage probabilities --------------------
# # following function is for creating coverage probability plots
quantile_df <- function(x, probs = c(0.25, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = c('lower', 'upper')
  ) |> pivot_wider(names_from = quant, values_from = val)
}

read_post_dep <- function(gauge, dep_level, lhood, data_num) {
  fit <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs", gauge, gauge, lhood, dep_level, data_num))$samples
  return(as.numeric(fit[,"dep"]))
}

create_tib_ci <- function(post_dep, ci_level, truth) {
  alpha <- (1-ci_level)/2
  probs_vec <- c(alpha, 1 - alpha)
  temp <- quantile_df(post_dep, probs = probs_vec)
  return(temp |> mutate(level = ci_level,
                        truth = truth,
                        coverage = ifelse(lower <= truth & upper >= truth, 1, 0)))
}

coverage_list <- function(gauge, dep_level, lhood, data_num, truth) {
  post_dep <- read_post_dep(gauge, dep_level, lhood, data_num)
  levels_vec <- seq(0.05, 0.95, by = 0.05)
  cov_list <- sapply(levels_vec, function(x) create_tib_ci(post_dep, x, truth), simplify = FALSE) |> 
    bind_rows() |>
    mutate(dataset = data_num)
  return(list(post_dep = post_dep, cov_list = cov_list))
}

extract_coverage <- function(gauge, lhood) {
  temp_list <- vector("list", 3)
  names(temp_list) <- c("high", "mid", "low")
  if (gauge == "gauss") {
    true_vals <- c(0.9, 0.5, 0.1)
    for(i in seq_along(true_vals)) {
      temp_list[[i]] <- sapply(1:200, function(y) coverage_list("gauss", names(temp_list)[i], lhood, y, true_vals[i]))
    }
  } else {
    true_vals <- c(0.1, 0.5, 0.9)
    for(i in seq_along(true_vals)) {
      temp_list[[i]] <- sapply(1:200, function(y) coverage_list("logistic", names(temp_list)[i], lhood, y, true_vals[i]))
    }
  }
  cov_plot_df <- lapply(temp_list, function(x) x[2,] |> bind_rows()) |>
    bind_rows() |>
    group_by(level, truth) |>
    summarize(p_hat = mean(coverage),
              sd = sd(coverage)) |>
    ungroup() |>
    mutate(se = sd/sqrt(200),
           lb = p_hat - qnorm(0.975) * se,
           ub = p_hat + qnorm(0.975) * se,
           truth = as.factor(truth))
  # med_mean_df <- lapply(temp_list, function(x) lapply(x[1,], function(y)  c("median" = median(y), "mean" = mean(y)))) |> 
  #   bind_rows() |> 
  #   mutate(truth = rep(true_vals, times = rep(200, 3)))
  # return(list(cov_plot_df = cov_plot_df,
  #             med_mean_df = med_mean_df))
  obj_name <- sprintf("%s_%s_coverage", gauge, lhood)
  assign(obj_name, cov_plot_df, envir = .GlobalEnv)
}

dep_lhood_combos <- expand_grid(dep_type = c("gauss", "logistic"),
                                lhood = c("trunc", "cens"))

with_progress({
  # Create a progress handler
  p <- progressor(steps = nrow(dep_lhood_combos))
  
  apply(dep_lhood_combos, 1, function(row) {
    extract_coverage(gauge = row["dep_type"],
                     lhood = row["lhood"])
    p()
  })
})

## Create boxplots and coverage plots ----------
plot_coverage <- function(cov_tibble, true_dep) {
  cov_tibble |>
    filter(truth == true_dep) |>
    ggplot(aes(x = level, y = p_hat, color = lhood)) +
    geom_abline(intercept = 0, slope = 1, color = "black") +
    geom_point(position = position_dodge(width = 0.01)) +
    geom_linerange(aes(ymin = lb, ymax = ub),
                   position = position_dodge(width = 0.01),
                   linetype = 2) +
    scale_color_grafify("r4", ColSeq = FALSE) +
    theme_classic() +
    ylim(-0.01, 1.01) +
    xlim(-0.01, 1.01) +
    xlab("Nominal Rate") +
    ylab("Empirical Rate") +
    labs(color = "Likelihood") +
    theme(panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = 'transparent'),
      axis.text = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.2)),
      legend.title = element_text(size = rel(1.2)),
      legend.text = element_text(size = rel(1.2)),
      legend.background = element_rect(fill = 'transparent', color = 'transparent'),
      legend.position = "inside",
      legend.position.inside = c(0.2, 0.85))
}

# plot_coverage(gauss_all_coverage, 0.9)

gauss_high_boxplot + plot_coverage(gauss_all_coverage, 0.9) + plot_layout(widths = c(0.95, 0.95))
ggsave("figures/gauss_high_coverage.pdf",
       bg = 'transparent',
       dpi = 320,
       width = 8,
       height = 4)
knitr::plot_crop("figures/gauss_high_coverage.pdf")

gauss_mid_boxplot + plot_coverage(gauss_all_coverage, 0.5) + plot_layout(widths = c(0.95, 0.95))
ggsave("figures/gauss_mid_coverage.pdf",
       bg = 'transparent',
       dpi = 320,
       width = 8,
       height = 4)
knitr::plot_crop("figures/gauss_mid_coverage.pdf")

gauss_low_boxplot + plot_coverage(gauss_all_coverage, 0.1) + plot_layout(widths = c(0.95, 0.95))
ggsave("figures/gauss_low_coverage.pdf",
       bg = 'transparent',
       dpi = 320,
       width = 8,
       height = 4)
knitr::plot_crop("figures/gauss_low_coverage.pdf")

logistic_high_boxplot + plot_coverage(logistic_all_coverage, 0.1) + plot_layout(widths = c(0.95, 0.95))
ggsave("figures/logistic_high_coverage.pdf",
       bg = 'transparent',
       dpi = 320,
       width = 8,
       height = 4)
knitr::plot_crop("figures/logistic_high_coverage.pdf")

logistic_mid_boxplot + plot_coverage(logistic_all_coverage, 0.5) + plot_layout(widths = c(0.95, 0.95))
ggsave("figures/logistic_mid_coverage.pdf",
       bg = 'transparent',
       dpi = 320,
       width = 8,
       height = 4)
knitr::plot_crop("figures/logistic_mid_coverage.pdf")

logistic_low_boxplot + plot_coverage(logistic_all_coverage, 0.9) + plot_layout(widths = c(0.95, 0.95))
ggsave("figures/logistic_low_coverage.pdf",
       bg = 'transparent',
       dpi = 320,
       width = 8,
       height = 4)
knitr::plot_crop("figures/logistic_low_coverage.pdf")

qsave(gauss_all_coverage, "figures/calibration_qs/gauss_coverage.qs")
qsave(logistic_all_coverage, "figures/calibration_qs/logistic_coverage.qs")




