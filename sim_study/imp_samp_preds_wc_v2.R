# =============================================================================
# Importance-sampling predictions using the Wadsworth-Campbell (wc) variable-k
# extrapolation scheme: k is chosen adaptively from the gauge threshold
# relative to each prediction box, so the importance samples overlap the box
# rather than being drawn at the observed scale. ("wc" throughout this script
# refers to Wadsworth-Campbell.) Supersedes imp_samp_preds_wc.R, which had a
# bug in gen_new_data().
#
# Run after all MCMC, loglik, weight extraction, and parameter extraction steps.
#
# Inputs:    fits_and_weights/post_params_joint/...qs
#            fits_and_weights/wts_joint_model/...qs
#            data/{dep_type}/{dep_level}_{i}.json
# Outputs:   figures/is_preds_boxplots/joint/{dep_type}_{dep_level}_{box}_sim_angles_v2.pdf
#            figures/is_preds_boxplots/joint/plot_objects/{...}_v2.qs
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

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

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

est_volume <- function(n = 100, pars, gauge_type) {
  temp <- seq(0, 1, length.out = n)
  grid <- as.matrix(expand.grid(temp, temp))
  gauge_fn <- get_gauge_function(gauge_type)
  gx <- gauge_fn(grid[,1], grid[,2], pars)
  return(mean(gx <= 1))
}

asym_log_gauge <- function(w1, w2, dep_par = 0.5) {
  r_inv <- 1/dep_par
  return(pmin((w1 + w2), (r_inv * pmax(w1, w2) + (1-r_inv)*pmin(w1,w2))))
}



temp <- seq(0, 1, length.out = 100)
grid <- as.matrix(expand.grid(temp, temp))
gx <- asym_log_gauge(grid[,1], grid[,2], 0.5)
mean(gx <= 1)

est_volume(n = 100, pars = 0.9, "asym_log")
w <- seq(0, 1, length.out = 300)
gw <- asym_log_gauge(w, 1-w, 0.9)
plot(w/gw, (1-w)/gw)
asym_log_gauge


star_dens <- function(w1, pars, gauge_type) {
  w2 <- 1 - w1
  mc_star <- est_volume(n = 100, pars, gauge_type)
  gauge_fn <- get_gauge_function(gauge_type)
  gw <- gauge_fn(w1, w2, pars)
  return(1 / (gw^2 * 2 * mc_star))
}

mix_dens_rng <- function(size = 1000, mean_params) {
  alphas <- as.numeric(mean_params[grepl("alphastar", names(mean_params))])
  betas <- as.numeric(mean_params[grepl("betastar", names(mean_params))])
  weights <- as.numeric(mean_params[grepl("probs", names(mean_params))])
  n <- length(weights)
  
  # randomly select which component to generate from, based on posterior probabilities
  dens_comp <- sample(1:n, size = size, TRUE, weights)
  angles <- rep(NA, size)
  for(i in 1:size) {
    angles[i] <- rbeta(1, shape1 = alphas[dens_comp[i]], shape2 = betas[dens_comp[i]])
  }
  return(angles)
}

imp_weights <- function(k, w_star, r0w_star, post_radial_alpha, post_rate) {
  numer <- pgamma(k * r0w_star, shape = post_radial_alpha, rate = post_rate, 
                  lower.tail = FALSE, log.p = TRUE)
  denom <- pgamma(r0w_star, shape = post_radial_alpha, rate = post_rate, 
                  lower.tail = FALSE, log.p = TRUE)
  return(exp(numer - denom))
}



gen_new_data <- function(k = 1, N = 50000, post_angular_params, post_radial_params, mix = TRUE, gauge_type) {
  # grab gauge function
  gauge_fn <- get_gauge_function(gauge_type)
  
  post_radial_dep <- as.numeric(post_radial_params[2:(length(post_radial_params) - 1)])
  post_radial_alpha <- post_radial_params[["alpha"]]
  
  if(mix) {
    sim_w <- mix_dens_rng(size = N, post_angular_params)
  } else {
    sim_w <- reject_sampling(n = N, post_radial_dep, gauge_type = gauge_type, dim = 2)
  }
  
  gw <- gauge_fn(sim_w, 1 - sim_w, post_radial_dep)
  r0w <- qgamma(0.95, post_radial_alpha, gw)
  
  if(k != 1) {
    wts <- imp_weights(k = k, w_star = sim_w, r0w_star = r0w, post_radial_alpha = post_radial_alpha, post_rate = gw)
    resample_idx <- sample(1:N, size = N, replace = TRUE, prob = wts)
    sim_w <- sim_w[resample_idx]
    gw <- gw[resample_idx]
    r0w <- k * r0w[resample_idx]
  }
  
  sim_r <- qgamma(1 - runif(N) * pgamma(r0w, shape = post_radial_alpha, 
                                        rate = gw, lower.tail = F),
                  shape = post_radial_alpha, rate = gw)
  xy <- cbind(sim_r * sim_w, sim_r * (1 - sim_w)) |> as_tibble() |> rename(x = V1, y = V2)
  
  return(list(xy = xy, r_over_k = ifelse(k != 1, mean(wts), 1)))
}

# probability calculation
is_prob_pred <- function(data_obs,
                         post_radial_params, 
                         post_angular_params,
                         box = "b1",
                         gauge_type, 
                         mix = TRUE,
                         k_var = FALSE) {
  # specify dimensions of box
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  # grab gauge function
  gauge_fn <- get_gauge_function(gauge_type)
  
  post_radial_dep <- as.numeric(post_radial_params[2:(length(post_radial_params) - 1)])
  post_radial_alpha <- post_radial_params[["alpha"]]
  
  if(k_var) {
    # determine threshold for the fitted gauge with posterior parameters
    gw_fitted <- gauge_fn(data_obs$W, 1 - data_obs$W, post_radial_dep)
    ctau_fitted <- quantile(gw_fitted * data_obs$R, 0.95)
    
    # create fake data to use in determining k value
    pseudo_pred <- expand_grid(x1_pseudo = seq(dim1[1], dim1[2], length.out=15), 
                               x2_pseudo = seq(dim2[1], dim2[2], length.out=15)) |> 
      mutate(r_pseudo = x1_pseudo + x2_pseudo,
             w1_pseudo = x1_pseudo / r_pseudo,
             w2_pseudo = x2_pseudo / r_pseudo)
    # determine ideal value of k using the above
    gw_pseudo <- gauge_fn(pseudo_pred$w1_pseudo, pseudo_pred$w2_pseudo, post_radial_dep)
    poss_k <- pseudo_pred$r_pseudo * gw_pseudo / ctau_fitted
    k <- max(round(min(poss_k), 1) - 0.5, 1)
  } else {
    k <- 1
  }
  
  sim_data <- gen_new_data(k = k, 
                           post_angular_params = post_angular_params, 
                           post_radial_params = post_radial_params, 
                           mix = mix, 
                           gauge_type = gauge_type)
  
  # subset IS samples based on box
  in_box <- which(
    with(
      sim_data$xy, 
      between(x, dim1[1], dim1[2]) & between(y, dim2[1], dim2[2]))
  )
  
  prob_pred <- length(in_box) / nrow(sim_data$xy) * sim_data$r_over_k * length(data_obs$idx)/data_obs$N
  
  return(list(prob_pred = prob_pred, in_box = in_box, k = k, xy = sim_data$xy))
}

preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, box, k_var) {
  post_radial <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_%s_radial.qs",
                               gauge, dep_type, dep_level, likelihood))
  post_ang_mix <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_ang_mix.qs",
                                dep_type, dep_level))
  post_ang_star <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_ang_star.qs",
                                 dep_type, dep_level, gauge))
  
  mix <- map_dbl(1:200, function(i) {
    is_prob_pred(data_obs = fload(sprintf("data/%s/%s_%s.json", dep_type, dep_level, i)),
                 post_radial_params = post_radial[i, ],
                 post_angular_params = post_ang_mix[i, ],
                 box = box,
                 gauge_type = gauge, 
                 mix = TRUE,
                 k_var = k_var)$prob_pred
  })
  
  star <- map_dbl(1:200, function(i) {
    is_prob_pred(data_obs = fload(sprintf("data/%s/%s_%s.json", dep_type, dep_level, i)),
                 post_radial_params = post_radial[i, ],
                 post_angular_params = post_ang_star[i, ],
                 box = box,
                 gauge_type = gauge, 
                 mix = FALSE,
                 k_var = k_var)$prob_pred
  })
  
  preds <- cbind(mix, star) |> as_tibble() |>
    mutate(method = gauge,
           dataset = 1:200)
  return(preds)
}

preds_by_dep_level_lhood <- function(dep_type, dep_level, likelihood, box, k_var) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(lapply(gauge_library, function(x) preds_by_gauge(gauge = x, 
                                                          dep_type = dep_type, 
                                                          dep_level = dep_level, 
                                                          likelihood = likelihood, 
                                                          box = box,
                                                          k_var = k_var)) |> 
           bind_rows())
}

weighted_preds_by_lhood <- function(dep_type, dep_level, likelihood, box, k_var) {
  preds <- preds_by_dep_level_lhood(dep_type = dep_type, 
                                    dep_level = dep_level, 
                                    likelihood = likelihood, 
                                    box = box, 
                                    k_var = k_var) |>
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

weighted_preds_by_level <- function(dep_type, dep_level, box, k_var) {
  lhood <- c("trunc", "cens")
  all_wts <- lapply(lhood, function(x) weighted_preds_by_lhood(dep_type = dep_type, 
                                                               dep_level = dep_level, 
                                                               likelihood = x,
                                                               box = box,
                                                               k_var = k_var)) |> bind_rows()
  return(all_wts)
}

# create function that wraps everything to gether to make boxplot of predictions
create_predictions_boxplot <- function(dep_type, dep_level, box, k_var) {
  
  dim1 <- c(10, 12)
  if(box == "b1") {
    dim2 <- dim1
  } else if(box == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  plot_filename <- sprintf("figures/is_preds_boxplots/joint/%s_%s_%s_sim_angles_v2.pdf",
                           dep_type, dep_level, box)
  qs_filename <-  sprintf("figures/is_preds_boxplots/joint/plot_objects/%s_%s_%s_sim_angles_v2.qs",
                          dep_type, dep_level, box)
  plot_title <- sprintf("%s, %s, (%s) x (%s)",
                        dep_type, dep_level, paste(dim1, collapse = ","), paste(dim2, collapse = ","))
  
  preds_tib <- weighted_preds_by_level(dep_type = dep_type, 
                                       dep_level = dep_level, 
                                       box = box,
                                       k_var = k_var)
  
  # True probability under the data-generating dependence parameters. Only the
  # gauss and logistic branches are exercised here (the Wadsworth-Campbell
  # comparison is run on those two dep_types); the Husler-Reiss branch is kept
  # for completeness with the same values used in data generation (0.1/1/3).
  true_prob <- get_true_prob(dep_type, dep_level, dim1, dim2)
  
  # create boxplot
  plot <- preds_tib |> ggplot(aes(x = ang_dens, y = preds, fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = true_prob, col = "darkgrey", linetype = "longdash") +
    theme_classic() +
    ggtitle(plot_title) + 
    # scale_fill_discrete(breaks = ~ .x[!is.na(.x)]) +
    xlab("Likelihood") + ylab("Prediction probabilities") + labs(fill = "")
  ggsave(plot_filename,
         plot = plot,
         bg = 'transparent',
         width = 8,
         height = 7,
         dpi = 320)
  qsave(plot, qs_filename)
  print(paste0(plot_filename, " has been saved"))
}

dep_types <- c("gauss", "logistic")
dep_levels <- c("high", "mid", "low")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)

with_progress({
  # Create a progress handler
  # p <- progressor(steps = nrow(all_combos))
  p <- progressor(steps = nrow(all_combos))
  
  # Apply the function using apply and update the progress bar
  apply(all_combos, 1, function(row) {
    create_predictions_boxplot(dep_type = row["dep_types"],
                               dep_level = row["dep_levels"],
                               box = row["boxes"],
                               k_var = TRUE)
    p()  # Update the progress bar
  })
})
