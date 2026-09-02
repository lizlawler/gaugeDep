# =============================================================================
# Full-posterior importance-sampling prediction on the real fire weather data:
# runs predictions across all posterior iterations (not just posterior medians)
# using parallelism via furrr. Produces the prediction uncertainty intervals
# shown in the final figures.
#
# Inputs:    fits_and_weights/post_params_joint/...qs (via extract_post_params_real_data.R)
#            fits_and_weights/wts_joint_model/...qs
#            data/raw/{data_type}_expo.qs
# Outputs:   real_data_preds/{data_type}_b{1,2,3}_all_iter.qs
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
source("extraction_scripts/extract_post_params_real_data.R")

library(furrr)
plan(multisession, workers = 4)
# 
options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

# data_type <- "redstone"

gauge_functions <- list(
  gauss = gauss_gauge,
  inv_log = inv_log_gauge,
  rectangular = rectangular_gauge,
  logistic = logistic_gauge,
  asym_log = asym_log_gauge,
  dirichlet = dirichlet_gauge
)

# Grab gauge function by string
get_gauge_function <- function(type_str) {
  if (!type_str %in% names(gauge_functions)) {
    stop("Unknown gauge type: ", type_str)
  }
  return(gauge_functions[[type_str]])
}

# Gamma density truncated below at xmin (the radial threshold), computed in
# log space then exponentiated.
trunc_gamma <- function(x, xmin, alpha, beta) {
  unnorm_pdf <- dgamma(x, shape = alpha, rate = beta, log = TRUE)
  norm_cst <- pgamma(xmin, shape = alpha, rate = beta, lower.tail = F, log.p = TRUE)
  return(exp(unnorm_pdf - norm_cst))
}

# Monte Carlo estimate of the unit star-body volume {x : g(x) <= 1}, used to
# normalise the star-shaped angular density.
est_volume <- function(n = 100, pars, gauge_type) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  gauge_fn <- get_gauge_function(gauge_type)
  gx <- gauge_fn(grid[,1], grid[,2], pars)
  return(mean(gx <= 1))
}

# Star-shaped angular density f(w) = 1 / (g(w)^2 * 2 * vol(L)).
star_dens <- function(w1, pars, gauge_type) {
  w2 <- 1 - w1
  mc_star <- est_volume(n = 100, pars, gauge_type)
  gauge_fn <- get_gauge_function(gauge_type)
  gw <- gauge_fn(w1, w2, pars)
  return(1 / (gw^2 * 2 * mc_star))
}

# Angular mixture density: weighted sum of Beta components (NIMBLE fit).
mix_dens <- function(w, mean_params) {
  alphas <- as.numeric(mean_params[grepl("alphastar", names(mean_params))])
  betas <- as.numeric(mean_params[grepl("betastar", names(mean_params))])
  weights <- as.numeric(mean_params[grepl("probs", names(mean_params))])
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}

# Importance samples: a bivariate Gaussian centred on the target box with
# covariance scaled to the box's width, restricted to the positive quadrant
# and converted to polar (r, w) form.
gen_is_samples <- function(dim1, dim2, total_n = 6000) {
  # generate samples from importance distribution
  is_samp_mvn <- mvtnorm::rmvnorm(total_n, 
                                  mean = c(mean(dim1), mean(dim2)), 
                                  sigma = matrix(c(diff(dim1), 0, 0, diff(dim2)), nrow = 2)) |>
    as_tibble() |>
    rename(x1 = V1, x2 = V2) |>
    filter(x1 >= 0, x2 >= 0) |> # adjust for truncation in later steps
    mutate(r = x1 + x2, w1 = x1 / r, w2 = x2 / r)
  
  return(is_samp_mvn)
}

# probability calculation
# Importance-sampling estimate of P(X in box) for one posterior draw. The
# joint (R, W) density is evaluated at the samples, weighted by
# target/proposal, and averaged over those landing in the box. The proposal
# density is renormalised by pmvnorm() to account for the positive-quadrant
# truncation, and the trailing 0.05 is the tail probability above the 95%
# threshold.
is_prob_pred <- function(imp_samples = NULL, 
                         post_radial_params, 
                         post_angular_params,
                         dim1,
                         dim2,
                         gauge_type, 
                         ang_dens = "star") {
  
  if(is.null(imp_samples)) {
    imp_samples <- gen_is_samples(dim1 = dim1, dim2 = dim2)
  }
  
  # grab gauge function
  gauge_fn <- get_gauge_function(gauge_type)
  
  post_radial_dep <- as.numeric(post_radial_params[2:length(post_radial_params)])
  post_radial_alpha <- post_radial_params[["alpha"]]
  
  # compute angular density
  ang_dens <- if (ang_dens == "star") {
    post_angular_dep <- as.numeric(post_angular_params[1:length(post_angular_params)])
    star_dens(imp_samples$w1, post_angular_dep, gauge_type)
  } else {
    mix_dens(imp_samples$w1, post_angular_params)
  }
  
  # estimate RW density
  gauge_vals <- gauge_fn(imp_samples$w1, imp_samples$w2, post_radial_dep)
  r0w <- qgamma(0.95, shape = post_radial_alpha, rate = gauge_vals)
  r_giv_w_dens <- trunc_gamma(imp_samples$r, r0w, alpha = post_radial_alpha, beta = gauge_vals)
  
  rw_dens <- r_giv_w_dens * ang_dens
  is_dens <- mvtnorm::dmvnorm(imp_samples[, 1:2], 
                              mean = c(mean(dim1), mean(dim2)), 
                              sigma = matrix(c(diff(dim1), 0, 0, diff(dim2)), nrow = 2)) * 
    imp_samples$r /
    mvtnorm::pmvnorm(lower = c(0,0),
                     mean = c(mean(dim1), mean(dim2)),
                     sigma = matrix(c(diff(dim1), 0, 0, diff(dim2)), nrow = 2)) # account for truncation above zero
  wts <- rw_dens/is_dens
  
  # indicator of being in box
  idx_in_box <- which(
    with(
      imp_samples,
      between(r, dim1[1] / w1, dim1[2] / w1) & 
        between(r, dim2[1] / (1-w1), dim2[2] / (1-w1)))
  )
  
  return(sum(wts[idx_in_box]) / nrow(imp_samples) * 0.05)
}

# Predictions for one gauge across every posterior iteration, under both the
# mixture and star angular densities. Box dimensions are station-specific and
# read from pred_boxes.qs (written by pred_task_real_data.R, which must run
# first).
preds_by_gauge <- function(gauge, likelihood, data, box) {
  post_radial <- extract_post_params_radial(gauge, likelihood, data, FALSE)
  post_ang_star <- extract_post_params_ang_star(gauge, data, FALSE)
  post_ang_mix <- extract_post_params_ang_mix(data, FALSE)
  
  box_dims <- qs::qread("data_analyses/pred_boxes.qs")[[data]][[box]]
  dim1 <- box_dims$dim1
  dim2 <- box_dims$dim2
  
  is_samp <- gen_is_samples(dim1 = dim1, dim2 = dim2)
  
  mix <- map_dbl(1:nrow(post_radial), function(i) {
    is_prob_pred(imp_samples = is_samp,
                 post_radial_params = post_radial[i, ],
                 post_angular_params = post_ang_mix[i, ],
                 dim1 = dim1, 
                 dim2 = dim2,
                 gauge_type = gauge, 
                 ang_dens = "mix")
  })
  
  star <- map_dbl(1:nrow(post_radial), function(i) {
    is_prob_pred(imp_samples = is_samp,
                 post_radial_params = post_radial[i, ],
                 post_angular_params = post_ang_star[i, ],
                 dim1 = dim1, 
                 dim2 = dim2,
                 gauge_type = gauge, 
                 ang_dens = "star")
  })
  
  preds <- cbind(mix, star) |> as_tibble() |>
    mutate(method = gauge, iter = 1:nrow(post_radial))
  return(preds)
}

preds_by_lhood <- function(likelihood, data, box, p) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  preds <- lapply(gauge_library, function(x) {
    temp <- preds_by_gauge(gauge = x,  
                           likelihood = likelihood, 
                           data = data,
                           box = box)
    p()
    temp
  }) |> bind_rows()
  return(preds)
}

# Combine the per-gauge predictions into BMA-weighted predictions (stacking
# and pseudo-BMA), per posterior iteration.
weighted_preds_by_lhood <- function(likelihood, data, box, p) {
  preds <- preds_by_lhood(likelihood = likelihood, 
                          data = data,
                          box = box,
                          p = p) |>
    pivot_longer(cols = c(mix, star), names_to = "ang_dens", values_to = "preds")
  wts_star <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_star.qs",
                            data, likelihood)) |> mutate(ang_dens = "star")
  wts_mix <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_mix.qs",
                           data, likelihood)) |> mutate(ang_dens = "mix")
  wts <- rbind(wts_star, wts_mix)
  wtd_preds <- suppressMessages(preds |> left_join(wts) |>
                                  mutate(stacking_preds = preds * stacking,
                                         pseudo_boot = pseudobma_boot * preds,
                                         pseudo_noboot = pseudobma_noboot * preds) |>
                                  group_by(iter, ang_dens) |>
                                  summarize(stacking_predictions = sum(stacking_preds),
                                            pseudobma_boot_preds = sum(pseudo_boot),
                                            pseudobma_noboot_preds = sum(pseudo_noboot)) |>
                                  ungroup())
  boxplot_wts <- wtd_preds |> 
    pivot_longer(cols = -c(iter, ang_dens), names_to = "method", values_to = "preds") |>
    mutate(method = case_when(grepl("stacking", method) ~ 'Stacking',
                              grepl("noboot", method) ~ 'Pseudo-BMA',
                              grepl("boot", method) ~ 'Pseudo-BMA+'),
           method = as.factor(method),
           ang_dens = paste0(likelihood, ", ", ang_dens))
  return(boxplot_wts)
}

weighted_preds <- function(data, box, p) {
  lhood <- c("trunc", "cens")
  all_wts <- lapply(lhood, function(x) weighted_preds_by_lhood(likelihood = x,
                                                               data = data,
                                                               box = box,
                                                               p = p)) |> 
    bind_rows()
  qsave(all_wts, sprintf("real_data_preds/%s_%s_all_iter.qs", data, box))
  # print(sprintf("Predictions for all iterations have been saved for box: %s", box))
}

data_type <- c("redstone", "friendmtn")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(data_type, boxes)

with_progress({
  p <- progressor(steps = nrow(all_combos) * 2 * 6)
  
  # Apply the function using apply and update the progress bar
  future_pmap(all_combos, function(data_type, boxes) {
    weighted_preds(data = data_type,
                   box = boxes,
                   p = p)
  }, .options = furrr_options(seed = TRUE))
})

