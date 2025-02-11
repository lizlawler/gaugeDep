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

# Function to get a gauge function by name
get_gauge_function <- function(type_str) {
  if (!type_str %in% names(gauge_functions)) {
    stop("Unknown gauge type: ", type_str)
  }
  return(gauge_functions[[type_str]])
}

trunc_gamma <- function(x, xmin, alpha, beta) {
  unnorm_pdf <- dgamma(x, shape = alpha, rate = beta, log = TRUE)
  norm_cst <- pgamma(xmin, shape = alpha, rate = beta, lower.tail = F, log.p = TRUE)
  return(exp(unnorm_pdf - norm_cst))
}

est_volume <- function(n = 100, pars, gauge_type) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  gauge_fn <- get_gauge_function(gauge_type)
  gx <- gauge_fn(grid[,1], grid[,2], pars)
  return(mean(gx <= 1))
}

star_dens <- function(w1, pars, gauge_type) {
  w2 <- 1 - w1
  mc_vol <- est_volume(n = 100, pars, gauge_type)
  gauge_fn <- get_gauge_function(gauge_type)
  gw <- gauge_fn(w1, w2, pars)
  return(1 / (gw^2 * 2 * mc_vol))
}

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

true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(mvtnorm::pmvnorm(lower = c(dim1_star[1],dim2_star[1]), upper = c(dim1_star[2],dim2_star[2]), corr = corr_matrix)[1])
}

true_bvevd_prob <- function(dim1, dim2, dep, model_type) {
  dim1_star <- qgev(pexp(dim1), loc = 0, scale = 1, shape = 0)
  dim2_star <- qgev(pexp(dim2), loc = 0, scale = 1, shape = 0)
  upper_right <- pbvevd(q = c(dim1_star[2], dim2_star[2]), model = model_type, dep = dep)
  upper_left <- pbvevd(q = c(dim1_star[1], dim2_star[2]), model = model_type, dep = dep)
  lower_right <- pbvevd(q = c(dim1_star[2], dim2_star[1]), model = model_type, dep = dep)
  lower_left <- pbvevd(q = c(dim1_star[1], dim2_star[1]), model = model_type, dep = dep)
  return(upper_right - upper_left - lower_right + lower_left)
}

gen_is_samples <- function(box = "b1", total_n = 5000) {
  # specify dimensions of box
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
    filter(x1 >= 0, x2 >= 0) |> # for B3, some samples may be < 0
    mutate(r = x1 + x2, w1 = x1 / r, w2 = x2 / r)
  
  return(is_samp_mvn)
}

# probability calculation
is_prob_pred <- function(imp_samples = NULL, 
                         post_radial_params, 
                         post_angular_params,
                         box = "b1",
                         renorm_gamma = TRUE, 
                         self_norm_is = FALSE,
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
  
  # subset IS samples based on box
  idx_in_box <- which(
    with(
      imp_samples, 
      r >= (10 / w1) & r <= (12 / w1) &
        r >= (dim2[1] / (1 - w1)) & r <= (dim2[2] / (1 - w1))
    )
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
  
  # estimate RW density
  gauge_vals <- gauge_fn(imp_samples$w1, imp_samples$w2, post_radial_dep)
  # r0w <- qgamma(0.95, shape = post_radial_params[["alpha"]], rate = gauge_vals)
  r_giv_w_dens <- if (renorm_gamma) {
    r0w <- qgamma(0.95, shape = post_radial_alpha, rate = gauge_vals)
    trunc_gamma(imp_samples$r, r0w, alpha = post_radial_alpha, beta = gauge_vals)
  } else {
    dgamma(imp_samples$r, shape = post_radial_alpha, rate = gauge_vals)
  }
  
  rw_dens <- r_giv_w_dens * ang_dens
  is_dens <- mvtnorm::dmvnorm(imp_samples[, 1:2], mean = c(mean(dim1), mean(dim2)), sigma = 2 * diag(2)) * imp_samples$r
  
  # numerator of vanilla IS and also self norm IS
  numer <- sum(rw_dens[idx_in_box] / is_dens[idx_in_box])
  # account for self normalization, if necessary
  denom <- if (self_norm_is) sum(rw_dens / is_dens) else nrow(imp_samples)
  # calculate fraction for IS part of probability calculation
  is_part <- numer / denom
  
  # account for renormalization of radii density
  return(is_part * ifelse(renorm_gamma, 0.05, 1))
}
# 
# test_radial <- qread("fits_and_weights/med_params_joint/dirichlet_gauss_high_trunc_radial.qs")
# test_ang_sb <- qread("fits_and_weights/med_params_joint/gauss_high_ang_sb.qs")
# test_ang_vol <- qread("fits_and_weights/med_params_joint/gauss_high_dirichlet_ang_vol.qs")
# 
# is_samps <- gen_is_samples("b1")
# is_prob_pred(post_radial_params = test_radial[1,], post_angular_params = test_ang_vol[1,], 
#              box = "b1", gauge_type = "dirichlet", ang_dens = "star")

preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, box) {
  post_radial <- qread(sprintf("fits_and_weights/med_params_joint/%s_%s_%s_%s_radial.qs",
                               gauge, dep_type, dep_level, likelihood))
  post_ang_sb <- qread(sprintf("fits_and_weights/med_params_joint/%s_%s_ang_sb.qs",
                               dep_type, dep_level))
  post_ang_star <- qread(sprintf("fits_and_weights/med_params_joint/%s_%s_%s_ang_vol.qs",
                                 dep_type, dep_level, gauge))
  
  mix <- map2_dbl(1:100, 1:100, 
                  ~ is_prob_pred(post_radial_params = post_radial[.x, ], 
                                 post_angular_params = post_ang_sb[.y, ],
                                 box = box, 
                                 gauge_type = gauge, ang_dens = "mix"))
  
  star <- map2_dbl(1:100, 1:100, 
                   ~ is_prob_pred(post_radial_params = post_radial[.x, ], 
                                  post_angular_params = post_ang_star[.y, ],
                                  box = box, 
                                  gauge_type = gauge, 
                                  ang_dens = "star"))
  preds <- cbind(mix, star) |> as_tibble() |>
    mutate(method = gauge,
           dataset = 1:100)
  # print(paste0("Predictions complete for gauge: ", gauge))
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
  # scenario <- sprintf("%s_%s_%s", 
  #                     dep_type, likelihood, dep_level)
  preds <- preds_by_dep_level_lhood(dep_type = dep_type, 
                                    dep_level = dep_level, 
                                    likelihood = likelihood, 
                                    box = box) |>
    pivot_longer(cols = c(mix, star), names_to = "ang_dens", values_to = "preds")
  wts_star <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_vol_%s.qs",
                            dep_type, likelihood, dep_level)) |> mutate(ang_dens = "star")
  wts_mix <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_sb_%s.qs",
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
  if(box == "b1") {
    dim2 <- dim1
  } else if(box == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  plot_filename <- sprintf("figures/is_preds_boxplots/joint/%s_%s_%s.pdf",
                           dep_type, dep_level, box)
  qs_filename <-  sprintf("figures/is_preds_boxplots/joint/plot_objects/%s_%s_%s.qs",
                          dep_type, dep_level, box)
  plot_title <- sprintf("%s, %s, (%s) x (%s)",
                        dep_type, dep_level, paste(dim1, collapse = ","), paste(dim2, collapse = ","))
  
  preds_tib <- weighted_preds_by_level(dep_type = dep_type, 
                                       dep_level = dep_level, 
                                       box = box)
  
  # determine true probability
  if(dep_type == "gauss") {
    levels_list <- list(low = 0.1, mid = 0.5, high = 0.9)
    true_prob <- true_gauss_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]))
  } else if(dep_type == "logistic"){
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "log")
  } else {
    levels_list <- list(low = 0.25, mid = 2, high = 6)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "hr")
  }
  
  # create boxplot
  plot <- preds_tib |> ggplot(aes(x = ang_dens, y = preds, fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = true_prob, col = "darkgrey", linetype = "longdash") +
    theme_classic() +
    ggtitle(plot_title) + 
    # scale_fill_discrete(breaks = ~ .x[!is.na(.x)]) +
    xlab("Likelihood and Threshold") + ylab("Prediction probabilities") + labs(fill = "")
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
  p <- progressor(steps = nrow(all_combos))
  
  # Apply the function using apply and update the progress bar
  apply(all_combos, 1, function(row) {
    p()  # Update the progress bar
    create_predictions_boxplot(dep_type = row["dep_types"],
                               dep_level = row["dep_levels"],
                               box = row["boxes"])
  })
})
