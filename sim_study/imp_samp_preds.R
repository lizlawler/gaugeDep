# =============================================================================
# Computes BMA-weighted importance-sampling predictions of bivariate exceedance
# probabilities for all 200 simulation study datasets under the level-set gauge
# model. Uses posterior-mean radial and angular parameters to evaluate the joint
# (R, W) density at importance samples drawn from a bivariate Gaussian centred
# on each prediction box, then applies stacking / pseudo-BMA weights across
# the 6 gauge functions and 2 angular density models (star, mix).
#
# Run after all MCMC, loglik calc, weight extraction, and parameter extraction
# steps have completed. Parallelises across scenarios via furrr.
#
# Inputs:    fits_and_weights/post_params_joint/...qs
#            fits_and_weights/wts_joint_model/...qs
# Outputs:   figures/is_preds_boxplots/joint/pred_tibbles/{dep_type}_{dep_level}_{box}_trunc.qs
#            figures/is_preds_boxplots/joint/mse_tables/{dep_type}_{dep_level}_{box}_trunc.qs
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
library(stringr)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

# Shared ground-truth probability functions
source("sim_study/true_prob_utils.R")

library(furrr)
plan(multisession, workers = parallel::detectCores()/2)

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

# Gamma density truncated below at xmin (the radial threshold), evaluated in
# log space for stability then exponentiated.
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

# Angular mixture density: weighted sum of Beta components (from NIMBLE fit).
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

# Draw importance samples from a bivariate Gaussian centred on the target
# box, keeping only the positive quadrant, and convert to polar (r, w) form.
gen_is_samples <- function(box = "b1", total_n = 5000) {
  # Box dimensions in Cartesian space
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  # generate samples from importance distribution
  is_samp_mvn <- mvtnorm::rmvnorm(total_n, c(mean(dim1), mean(dim2)), sigma = 2 * diag(2)) |>
    as_tibble() |>
    rename(x1 = V1, x2 = V2) |>
    filter(x1 >= 0, x2 >= 0) |> # adjust for this "truncation" in the is_prob_pred function
    mutate(r = x1 + x2, w1 = x1 / r, w2 = x2 / r)
  
  return(is_samp_mvn)
}

# Importance-sampling estimate of P(X in box) for one dataset and model.
# Evaluates the joint (R, W) density at the importance samples, weights each
# by (target density / proposal density), and averages over the samples that
# land inside the box. The trailing * 0.05 accounts for the tail region
# above the 95%% threshold.
is_prob_pred <- function(imp_samples = NULL, 
                         post_radial_params, 
                         post_angular_params,
                         box = "b1",
                         gauge_type, 
                         ang_dens = "star") {
  
  if(is.null(imp_samples)) {
    imp_samples <- gen_is_samples(box = box)
  }
  
  # grab gauge function
  gauge_fn <- get_gauge_function(gauge_type)
  
  # specify dimensions of box
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  post_radial_dep <- as.numeric(post_radial_params[2:(length(post_radial_params) - 1)])
  post_radial_alpha <- post_radial_params[["alpha"]]
  
  # compute angular density
  ang_dens <- if (ang_dens == "star") {
    post_angular_dep <- as.numeric(post_angular_params[1:(length(post_angular_params) - 1)])
    star_dens(imp_samples$w1, post_angular_dep, gauge_type)
  } else {
    mix_dens(imp_samples$w1, post_angular_params)
  }
  
  # Joint (R, W) density = truncated-Gamma radial density x angular density.
  # r0w is the 95%% radial quantile at each angle (the gauge threshold).
  gauge_vals <- gauge_fn(imp_samples$w1, imp_samples$w2, post_radial_dep)
  r0w <- qgamma(0.95, shape = post_radial_alpha, rate = gauge_vals)
  r_giv_w_dens <- trunc_gamma(imp_samples$r, r0w, alpha = post_radial_alpha, beta = gauge_vals)
  
  rw_dens <- r_giv_w_dens * ang_dens
  is_dens <- mvtnorm::dmvnorm(imp_samples[, 1:2], mean = c(mean(dim1), mean(dim2)), sigma = 2 * diag(2)) * 
    imp_samples$r / 
    mvtnorm::pmvnorm(lower = c(0,0), mean = c(mean(dim1), mean(dim2)), sigma = 2 * diag(2))[1]
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

preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, box) {
  post_radial <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_%s_radial.qs",
                               gauge, dep_type, dep_level, likelihood))
  post_ang_mix <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_ang_mix.qs",
                                dep_type, dep_level))
  post_ang_star <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_ang_star.qs",
                                 dep_type, dep_level, gauge))
  
  mix <- map_dbl(1:200, function(i) {
    is_prob_pred(post_radial_params = post_radial[i, ],
                 post_angular_params = post_ang_mix[i, ],
                 box = box,
                 gauge_type = gauge, 
                 ang_dens = "mix")
  })
  
  star <- map_dbl(1:200, function(i) {
    is_prob_pred(post_radial_params = post_radial[i, ],
                 post_angular_params = post_ang_star[i, ],
                 box = box,
                 gauge_type = gauge,
                 ang_dens = "star")
  })
  
  preds <- cbind(mix, star) |> as_tibble() |>
    mutate(method = gauge,
           dataset = 1:200)
  return(preds)
}

preds_by_dep_level_lhood <- function(dep_type, dep_level, likelihood, box) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(lapply(gauge_library, function(x) preds_by_gauge(gauge = x, 
                                                          dep_type = dep_type, 
                                                          dep_level = dep_level, 
                                                          likelihood = likelihood, 
                                                          box = box)) |> 
           bind_rows())
}

weighted_preds_by_lhood <- function(dep_type, dep_level, likelihood, box) {
  
  preds <- preds_by_dep_level_lhood(dep_type = dep_type, 
                                    dep_level = dep_level, 
                                    likelihood = likelihood, 
                                    box = box) |>
    pivot_longer(cols = c(mix, star), names_to = "ang_dens", values_to = "preds")
  wts_star <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_star_%s.qs",
                            dep_type, likelihood, dep_level)) |> mutate(ang_dens = "star")
  wts_mix <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_mix_%s.qs",
                           dep_type, likelihood, dep_level)) |> mutate(ang_dens = "mix")
  wts <- rbind(wts_star, wts_mix)
  wtd_preds <- suppressMessages(preds |> left_join(wts) |>
                                  mutate(stacking_preds = preds * stacking,
                                         pseudo_boot = pseudobma_boot * preds,
                                         pseudo_noboot = pseudobma_noboot * preds) |>
                                  group_by(dataset, ang_dens) |>
                                  summarize(stacking_predictions = sum(stacking_preds),
                                            pseudobma_boot_preds = sum(pseudo_boot),
                                            pseudobma_noboot_preds = sum(pseudo_noboot)) |>
                                  ungroup())
  boxplot_wts <- wtd_preds |> 
    pivot_longer(cols = -c(dataset, ang_dens), names_to = "method", values_to = "preds") |>
    mutate(method = case_when(grepl("stacking", method) ~ 'Stacking',
                              grepl("noboot", method) ~ 'Pseudo-BMA',
                              grepl("boot", method) ~ 'Pseudo-BMA+'),
           method = as.factor(method),
           ang_dens = paste0(likelihood, ", ", ang_dens))
  return(boxplot_wts)
}

weighted_preds_by_level <- function(dep_type, dep_level, box) {
  lhood <- c("trunc", "cens")
  all_wts <- lapply(lhood, function(x) weighted_preds_by_lhood(dep_type = dep_type, 
                                                               dep_level = dep_level, 
                                                               likelihood = x,
                                                               box = box)) |> bind_rows()
  return(all_wts)
}

# create function that wraps everything to gether to make boxplot of predictions
create_predictions_boxplot <- function(dep_type, dep_level, box) {
  
  dim1 <- c(10, 12)
  dim2 <- switch(box, b1 = dim1, b2 = c(6, 8), c(2, 4))

  preds_tib <- weighted_preds_by_level(dep_type = dep_type,
                                       dep_level = dep_level,
                                       box = box)
  qsave(preds_tib, sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_trunc.qs",dep_type, dep_level, box))
  
  true_prob <- get_true_prob(dep_type, dep_level, dim1, dim2)
  
  # create boxplot
  plot <- preds_tib |> ggplot(aes(x = ang_dens, y = pmax(preds, .Machine$double.eps), fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = pmax(true_prob, .Machine$double.eps), col = "darkgrey", linetype = "longdash") +
    scale_y_log10() +
    ggtitle(sprintf("%s, %s, (%s) x (%s)", dep_type, dep_level, paste(dim1, collapse = ","), paste(dim2, collapse = ","))) +
    theme_classic() +
    theme(panel.background = element_rect(fill='transparent', color='transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'),
          axis.text = element_text(size = rel(1.2)),
          axis.title = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1.2)),
          legend.title = element_text(size = rel(1.2)),
          legend.background = element_rect(fill='transparent', color='transparent')) +
    scale_fill_grafify(breaks = ~ .x[!is.na(.x)], palette = "r4", ColSeq = FALSE) + # use colorblind-friendly palette
    labs(fill = "BMA method") +
    scale_x_discrete(labels=c("Cens., Mix", "Cens., Star", "Trunc., Mix", "Trunc., Star")) +
    xlab("Likelihood, Angular Density") + ylab("Prediction probabilities")
  ggsave(sprintf("figures/is_preds_boxplots/joint/%s_%s_%s_trunc.pdf", dep_type, dep_level, box),
         plot = plot,
         bg = 'transparent',
         width = 8,
         height = 7,
         dpi = 320)
  qsave(plot, sprintf("figures/is_preds_boxplots/joint/plot_objects/%s_%s_%s_trunc.qs",dep_type, dep_level, box))
  print(sprintf("Boxplot for %s %s, in box %s has been saved", dep_level, dep_type, box))

  mse_table <- preds_tib |> 
    mutate(truth = true_prob,
           diff = preds - truth) |>
    group_by(method, ang_dens) |>
    summarise(mse = mean(diff^2)) |>
    ungroup() |>
    mutate(log_mse = log(mse), rmse_norm = sqrt(mse)/true_prob) |>
    arrange(rmse_norm)
  qsave(mse_table, sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_trunc.qs",dep_type, dep_level, box))
}

dep_types <- c("gauss", "logistic", "husler_reiss")
dep_levels <- c("high", "mid", "low")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)

with_progress({
  p <- progressor(steps = nrow(all_combos))
  
  future_pmap(all_combos, function(dep_types, dep_levels, boxes) {
    create_predictions_boxplot(dep_type = dep_types,
                               dep_level = dep_levels,
                               box = boxes)
    p()  # Update the progress bar
  }, .options = furrr_options(seed = TRUE))
})

