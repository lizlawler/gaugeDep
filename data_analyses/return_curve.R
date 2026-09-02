# =============================================================================
# Computes and plots return-level curves (joint exceedance probability as a
# function of return period) for the real fire weather data. Uses full-posterior
# importance-sampling predictions loaded from real_data_preds/.
#
# Inputs:    real_data_preds/{data_type}_b{1,2,3}_all_iter.qs
# Outputs:   figures/return_diag_plot_*.pdf
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(evd)
library(progressr)
library(RcppSimdJson)
library(gaugeDependence)
library(qs)
library(grafify)
library(patchwork)

## Uncomment these lines if returns have not been calculated yet ------
# source("extraction_scripts/extract_post_params_real_data.R")
# #
# options(rlib_name_repair_verbosity = "quiet")
# handlers("cli")
#
# library(furrr)
# # plan(multisession, workers = parallel::detectCores()/2)
# plan(multisession, workers = 6)
#
# gauge_functions <- list(
#   gauss = gauss_gauge,
#   inv_log = inv_log_gauge,
#   rectangular = rectangular_gauge,
#   logistic = logistic_gauge,
#   asym_log = asym_log_gauge,
#   dirichlet = dirichlet_gauge
# )
#
# # Grab gauge function by string
# get_gauge_function <- function(type_str) {
#   if (!type_str %in% names(gauge_functions)) {
#     stop("Unknown gauge type: ", type_str)
#   }
#   return(gauge_functions[[type_str]])
# }
#
# return_by_gauge_lhood_tau <- function(r, w, tau, dataname, gauge, likelihood) {
#   return_level <- 1 - 1/tau
#   gauge_fn <- get_gauge_function(gauge)
#
#   # read in posterior params
#   post_radial <- extract_post_params_radial(gauge = gauge, likelihood = likelihood, data = dataname, summarize = FALSE)
#   post_shape <- post_radial$alpha
#   post_dep <- post_radial[,-1]
#
#   # calculate posterior rate parameter
#   gauge_vals <- apply(post_dep, 1,function(row) gauge_fn(w, 1-w, as.numeric(row)))
#
#   # return level quantiles
#   r0w_post <- matrix(qgamma(return_level, shape = rep(post_shape, each = length(w)),
#                             rate = as.vector(gauge_vals)),
#                      nrow = length(w))
#   exceed_mat <- sweep(r0w_post, 1, r, "<")  # TRUE where r > r0w_post
#   est_returns <- colMeans(exceed_mat)
#
#   est_returns_tib <- tibble(
#     est_val = est_returns,
#     tau = tau,
#     iter = seq_along(est_returns),
#     likelihood = likelihood,
#     gauge = gauge
#   )
#   qsave(x = est_returns_tib, file = sprintf("figures/return_sets/%s_%s_%s_%s.qs", dataname, gauge, likelihood, tau))
#   rm(post_radial, gauge_vals, r0w_post, est_returns_tib)
#   gc()
# }
#
# data_type <- c("redstone", "friendmtn")
# gauge_library <- c("gauss", "inv_log", "rectangular", "logistic", "dirichlet", "asym_log")
# lhood <- c("cens", "trunc")
# tau_vals <- seq(10, 1000, by = 10)
# all_combos <- expand_grid(data_type, gauge_library, lhood, tau_vals)
#
# with_progress({
#   p <- progressor(steps = nrow(all_combos))
#
#   future_pmap(all_combos, function(data_type, gauge_library, lhood, tau_vals) {
#     data <- qread(sprintf("data/raw/%s_expo.qs", data_type))
#     return_by_gauge_lhood_tau(
#       r = data$R,
#       w = data$W,
#       tau = tau_vals,
#       dataname = data_type,
#       gauge = gauge_library,
#       likelihood = lhood
#     )
#     p()
#   }, .options = furrr_options(seed = TRUE))
# })
#
# full_est_return <- tibble()
# # data_type <- "friendmtn"
# all_files <- list.files(path = "figures/return_sets/", pattern = data_type, full.names = TRUE)
# for(i in seq_along(all_files)) {
#   full_est_return <- rbind(full_est_return, qread(all_files[i]))
# }
# qsave(full_est_return, sprintf("figures/return_sets/%s_all_est_returns.qs", data_type))

## start from here after saving objects above -------
create_full_return_set <- function(data_type, ang_dens) {
  full_est_return <- qread(sprintf("figures/return_sets/%s_all_est_returns.qs", data_type))
  wts_cens <- qread(sprintf("fits_and_weights/wts_joint_model/%s_cens_%s.qs", data_type, ang_dens)) |>
    rename(gauge = method) |>
    mutate(likelihood = "cens")
  wts_trunc <- qread(sprintf("fits_and_weights/wts_joint_model/%s_trunc_%s.qs", data_type, ang_dens)) |>
    rename(gauge = method) |>
    mutate(likelihood = "trunc")
  both_wts <- rbind(wts_cens, wts_trunc)

  full_est_return <- full_est_return |>
    left_join(both_wts) |>
    mutate(
      stacking_ests = est_val * stacking,
      pseudo_boot = pseudobma_boot * est_val,
      pseudo_noboot = pseudobma_noboot * est_val
    ) |>
    group_by(likelihood, tau, iter) |>
    summarize(
      stacking_est_return = sum(stacking_ests),
      pseudo_boot_est_return = sum(pseudo_boot),
      pseudo_noboot_est_return = sum(pseudo_noboot)
    ) |>
    ungroup()

  return(full_est_return |>
    pivot_longer(cols = -c(likelihood, tau, iter), names_to = "method", values_to = "ests") |>
    mutate(
      method = case_when(
        grepl("stacking", method) ~ "Stacking",
        grepl("noboot", method) ~ "Pseudo-BMA",
        grepl("boot", method) ~ "Pseudo-BMA+"
      ),
      method = as.factor(method),
      likelihood = as.factor(likelihood)
    ) |>
    group_by(likelihood, tau, method) |>
    summarize(
      median = median(ests),
      mean = mean(ests),
      lb = quantile(ests, 0.025),
      ub = quantile(ests, 0.975)
    ) |>
    ungroup() |>
    mutate(dataset = data_type))
}

friendmtn_returns_mix <- create_full_return_set("friendmtn", "mix")
friendmtn_returns_star <- create_full_return_set("friendmtn", "star")
redstone_returns_mix <- create_full_return_set("redstone", "mix")
redstone_returns_star <- create_full_return_set("redstone", "star")
full_est_return_mix <- rbind(friendmtn_returns_mix, redstone_returns_mix) |>
  mutate(dataset = case_when(
    dataset == "redstone" ~ "Redstone",
    dataset == "friendmtn" ~ "Friend Mountain"
  ))
full_est_return_star <- rbind(friendmtn_returns_star, redstone_returns_star) |>
  mutate(dataset = case_when(
    dataset == "redstone" ~ "Redstone",
    dataset == "friendmtn" ~ "Friend Mountain"
  ))

return_plot <- function(full_tib, stat) {
  full_tib |>
    filter(method == "Pseudo-BMA+") |>
    ggplot(aes(x = log(tau), y = -log(.data[[stat]]), color = likelihood, group = likelihood)) +
    geom_point(alpha = 0.8, size = 0.95) +
    scale_x_continuous(limits = c(2, 8), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(limits = c(2, 8), expand = expansion(mult = c(0, 0.01))) +
    geom_abline(slope = 1, intercept = 0, linetype = "solid") +
    geom_line(aes(x = log(tau), y = -log(lb), color = likelihood), linetype = "dashed") +
    geom_line(aes(x = log(tau), y = -log(ub), color = likelihood), linetype = "dashed") +
    scale_color_grafify(palette = "r4", labels = c("Censored", "Truncated"), ColSeq = FALSE) +
    facet_wrap(. ~ dataset) +
    ylab(expression("log(T"["est"] ~ ")")) +
    xlab(expression("log(T"["emp"] ~ ")")) +
    labs(color = "Likelihood") +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "transparent", color = "transparent"),
      plot.background = element_rect(fill = "transparent", color = "transparent"),
      axis.text = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.2)),
      legend.text = element_text(size = rel(1.2)),
      legend.title = element_text(size = rel(1.2)),
      legend.position = "inside",
      legend.position.inside = c(0.85, 0.2),
      panel.spacing.x = unit(1, "cm", data = NULL),
      panel.spacing.y = unit(0.5, "cm", data = NULL),
      # strip.text.x = element_text(size = rel(1.2)),
      strip.text = element_blank(),
      legend.background = element_rect(fill = "transparent", color = "transparent")
    )
}

## angular density = mixture
return_plot(full_est_return_mix, "median")
ggsave(
  filename = "figures/return_diag_plot_mix_median.pdf",
  dpi = 320,
  width = 10.3,
  height = 5,
  bg = "transparent"
)
knitr::plot_crop("figures/return_diag_plot_mix_median.pdf")

return_plot(full_est_return_mix, "mean")
ggsave(
  filename = "figures/return_diag_plot_mix_mean.pdf",
  dpi = 320,
  width = 10.3,
  height = 5,
  bg = "transparent"
)
knitr::plot_crop("figures/return_diag_plot_mix_mean.pdf")

## angular density = star-shaped
return_plot(full_est_return_star, "median")
ggsave(
  filename = "figures/return_diag_plot_star_median.pdf",
  dpi = 320,
  width = 10.3,
  height = 5,
  bg = "transparent"
)
knitr::plot_crop("figures/return_diag_plot_star_median.pdf")

return_plot(full_est_return_star, "mean")
ggsave(
  filename = "figures/return_diag_plot_star_mean.pdf",
  dpi = 320,
  width = 10.3,
  height = 5,
  bg = "transparent"
)
knitr::plot_crop("figures/return_diag_plot_star_mean.pdf")
