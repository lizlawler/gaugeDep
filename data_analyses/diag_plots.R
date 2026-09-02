# =============================================================================
# Produces PP/QQ diagnostic plots comparing the fitted marginal EGPD model
# to the observed fire weather data (ERC, FWI) for each station and gauge
# function. Loads posterior parameters via extract_post_params_real_data.R.
#
# Inputs:    fits_and_weights/post_params_joint/...qs (via extract_post_params_real_data.R)
#            data/raw/{data_type}_expo.qs
# Outputs:   figures/diag_plots/{data_type}_{gauge}_{likelihood}.qs
#            figures/pp_qq_diag_plots.pdf (and variants)
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
source("extraction_scripts/extract_post_params_real_data.R")


## Uncomment these lines if cdf by iteration for each dataset hasn't been calculated ----
# library(furrr)
# plan(multisession, workers = 3)
# 
# options(rlib_name_repair_verbosity = "quiet")
# handlers("cli")
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
# cdf_by_gauge <- function(data_type, gauge, likelihood) {
#   data <- qs::qread(sprintf("data/raw/%s_expo.qs", data_type))
#   idx <- data$idx
#   r <- data$R[idx]
#   w <- data$W[idx]
#   r0w <- data$r0_w[idx]
#   n0 <- as.numeric(data$n0)
#   
#   # read in posterior params
#   post_radial_params <- extract_post_params_radial(gauge, likelihood, data_type, FALSE)
#   post_dep <- post_radial_params[,-1]
#   post_shape <- post_radial_params$alpha
#   
#   # calculate rate parameter
#   gauge_fn <- get_gauge_function(gauge)
#   gauge_vals <- apply(post_dep, 1, function(row) gauge_fn(w, 1 - w, as.numeric(row)))
#   
#   # compute truncated cdf
#   num <- pgamma(r, shape = rep(post_shape, each = length(w)), rate = as.vector(gauge_vals), lower.tail = FALSE)
#   denom <- pgamma(r0w, shape = rep(post_shape, each = length(w)), rate = as.vector(gauge_vals), lower.tail = FALSE)
#   cdf <- 1 - (num / denom)
#   
#   cdf_tib <- cdf |> 
#     as_tibble() |> 
#     mutate(method = gauge, 
#            lhood = likelihood,
#            id = rep(1:n0, times = 4000), # every n0 rows is 1 iteration
#            iter = rep(1:4000, each = n0))
#   qsave(cdf_tib, sprintf("figures/diag_plots/%s_%s_%s.qs", data_type, gauge, likelihood))
#   rm(post_radial_params, post_shape, post_dep, gauge_vals, cdf, cdf_tib)
#   gc()
# }
# 
# data_names <- c("redstone", "friendmtn")
# gauge_library <- c("gauss", "inv_log", "rectangular", "logistic", "dirichlet", "asym_log")
# lhood <- c("cens", "trunc")
# all_combos <- expand_grid(data_names, gauge_library, lhood)
# 
# with_progress({
#   p <- progressor(steps = nrow(all_combos))
#   
#   future_pmap(all_combos, function(data_names, gauge_library, lhood) {
#     cdf_by_gauge(
#       data_type = data_names,
#       gauge = gauge_library,
#       likelihood = lhood
#     )
#     p()
#   }, .options = furrr_options(seed = TRUE))
# })
# 
# data_types <- c("redstone", "friendmtn")
# invisible(sapply(data_types, function(x) {
#   full_est_cdf <- tibble()
#   all_files <- list.files(path = "figures/diag_plots/", pattern = x, full.names = TRUE)
#   for(i in seq_along(all_files)) {
#     full_est_cdf <- rbind(full_est_cdf, qread(all_files[i]))
#   }
#   qsave(full_est_cdf, sprintf("figures/diag_plots/%s_all_est_cdf.qs", x))
# }))


## Create diagnostic plots -----
weighted_cdf_by_lhood <- function(data_type, ang_dens) {
  full_cdf_tib <- qread(sprintf("figures/diag_plots/%s_all_est_cdf.qs", data_type))
  wts_cens <- qread(sprintf("fits_and_weights/wts_joint_model/%s_cens_%s.qs", data_type, ang_dens)) |>
    mutate(lhood = "cens")
  wts_trunc <- qread(sprintf("fits_and_weights/wts_joint_model/%s_trunc_%s.qs", data_type, ang_dens)) |>
    mutate(lhood = "trunc")
  all_wts <- rbind(wts_cens, wts_trunc)
  wtd_cdf <- suppressMessages(full_cdf_tib |> left_join(all_wts) |>
                                mutate(stacking_preds = value * stacking,
                                       pseudo_boot = pseudobma_boot * value,
                                       pseudo_noboot = pseudobma_noboot * value) |>
                                group_by(id, lhood, iter) |>
                                summarize(stacking_predictions = sum(stacking_preds),
                                          pseudobma_boot_preds = sum(pseudo_boot),
                                          pseudobma_noboot_preds = sum(pseudo_noboot),
                                          .groups = "drop")) |>
    pivot_longer(cols = -c(id, lhood, iter), names_to = "method", values_to = "cdf") |>
    mutate(method = case_when(grepl("stacking", method) ~ "Stacking",
                              grepl("noboot", method) ~ "Pseudo-BMA",
                              grepl("boot", method) ~ "Pseudo-BMA+",
                              .default = method)) |>
    group_by(lhood, id, method) |>
    summarize(median = median(cdf),
              mean = mean(cdf),
              lb = quantile(cdf, 0.025),
              ub = quantile(cdf, 0.975),
              .groups = "drop") |>
    group_by(lhood, method) |>
    arrange(median, .by_group = TRUE) |>
    mutate(empir_prob = row_number() / (n() + 1),
           ulb = qbeta(0.025, row_number(), n() + 1 - row_number()),
           uub = qbeta(0.975, row_number(), n() + 1 - row_number())) |> 
    ungroup() |>
    mutate(dataset = data_type,
           lhood = case_when(lhood == "trunc" ~ "Truncated",
                             lhood == "cens" ~ "Censored",
                             .default = lhood))
  return(wtd_cdf)
}

data_names <- c("redstone", "friendmtn")
angs <- c("mix", "star")
data_ang_combos <- expand_grid(data_names, angs)
invisible(apply(data_ang_combos, 1, function(row) {
  temp <- weighted_cdf_by_lhood(row["data_names"], row["angs"])
  assign(sprintf("%s_%s", row["data_names"], row["angs"]), temp, envir = .GlobalEnv)
  rm(temp)
  gc()
}))

plot_by_margin <- function(full_tib, expo = FALSE, expo_limits = NULL) {
  if(expo) {
    full_tib <- full_tib |> mutate(across(where(is.double), qexp))
    limits <- if (is.null(expo_limits)) c(NA, 11) else expo_limits
  } else {
    limits <- c(NA, 1)
  }
  
  p <- full_tib |> filter(method == "Pseudo-BMA+") |>
    ggplot(aes(x = empir_prob, y = median, group = lhood, color = lhood)) +
    geom_point(alpha = 0.8, size = 1.2) + 
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_grafify(palette = "r4", ColSeq = FALSE) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = lhood), alpha = 0.35, color = NA) +
    scale_fill_grafify(palette = "r4", ColSeq = FALSE, guide = "none") +
    geom_line(aes(x = empir_prob, y = ulb), linetype = "dashed", color = "darkgrey") + 
    geom_line(aes(x = empir_prob, y = uub), linetype = "dashed", color = "darkgrey") + 
    labs(color = "Likelihood",
         y = "Model CDF",
         x = "Empirical CDF") +
    scale_x_continuous(limits = limits, expand = expansion(mult = c(0, 0.03))) +
    scale_y_continuous(limits = limits, expand = expansion(mult = c(0, 0.03))) +
    theme_classic() +
    theme(panel.background = element_rect(fill = 'transparent', color = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = 'transparent'),
          axis.text = element_text(size = rel(1.3)),
          axis.title = element_text(size = rel(1.3)))
  if(expo) {
    return(p + theme(legend.text = element_text(size = rel(1.3)),
                     legend.title = element_text(size = rel(1.3)),
                     legend.position = "inside",
                     legend.position.inside = c(0.8,0.2),
                     legend.background = element_rect(fill = 'transparent', color = 'transparent')) +
             guides(color = guide_legend(override.aes = list(size = 2, linewidth = 1))))
  } else {
    return(p + theme(legend.position = "none"))
  }
}

plot_by_data <- function(full_tib, expo_limits = NULL) {
  p <- plot_by_margin(full_tib, FALSE)
  q <- plot_by_margin(full_tib, TRUE, expo_limits)
  return(p + plot_spacer() + q + plot_layout(widths = c(0.95, 0.05, 0.95)))
}

## angular density = mixture
plot_by_data(friendmtn_mix, expo_limit = c(NA, 11)) / 
  plot_spacer() / 
  plot_by_data(redstone_mix, expo_limit = c(NA, 11.5)) + 
  plot_layout(heights = c(0.95, 0.05, 0.95))

ggsave("figures/pp_qq_diag_plots_mix.pdf",
       bg = 'transparent',
       width = 9,
       height = 9,
       dpi = 320)
knitr::plot_crop("figures/pp_qq_diag_plots_mix.pdf")

## angular density = star
plot_by_data(friendmtn_star, expo_limit = c(NA, 11)) / 
  plot_spacer() / 
  plot_by_data(redstone_star, expo_limit = c(NA, 11.5)) + 
  plot_layout(heights = c(0.95, 0.05, 0.95))

ggsave("figures/pp_qq_diag_plots_star.pdf",
       bg = 'transparent',
       width = 9,
       height = 9,
       dpi = 320)
knitr::plot_crop("figures/pp_qq_diag_plots_star.pdf")
