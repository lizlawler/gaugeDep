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

true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(mvtnorm::pmvnorm(lower = c(dim1_star[1],dim2_star[1]), 
                          upper = c(dim1_star[2],dim2_star[2]), 
                          corr = corr_matrix)[1])
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

imp_weights <- function(k, w, post_radial_dep, post_radial_alpha, gauge_type) {
  gauge_fn <- get_gauge_function(gauge_type)
  post_rate <- gauge_fn(w, 1 - w, post_radial_dep)
  r0w <- qgamma(0.95, post_radial_alpha, post_rate)
  numer <- pgamma(k * r0w, shape = post_radial_alpha, rate = post_rate, 
                  lower.tail = FALSE, log.p = TRUE)
  denom <- pgamma(r0w, shape = post_radial_alpha, rate = post_rate, 
                  lower.tail = FALSE, log.p = TRUE)
  return(list(wts = exp(numer - denom), r0w = r0w, rate = post_rate))
}

test_radial <- qread("fits_and_weights/post_params_joint/gauss_gauss_high_trunc_radial.qs")
test_ang_mix <- qread("fits_and_weights/post_params_joint/gauss_high_ang_mix.qs")
test_ang_star <- qread("fits_and_weights/post_params_joint/gauss_high_gauss_ang_star.qs")

new_w <- mix_dens_rng(50000, test_ang_mix[2,])
data <- RcppSimdJson::fload("data/gauss/high_2.json")
r <- data$R
w <- data$W
x <- r * w
y <- r * (1-w)
gw_fitted <- gauss_gauge(w, 1-w, test_radial$dep[2])
ctau_fitted <- quantile(gw_fitted * r, 0.95)
dim1 <- dim2 <- c(10,12)
# create fake data to use in determining k value
pseudo_pred <- expand_grid(x1_pseudo = seq(dim1[1], dim1[2], length.out=15), 
                           x2_pseudo = seq(dim2[1], dim2[2], length.out=15)) |> 
  mutate(r_pseudo = x1_pseudo + x2_pseudo,
         w1_pseudo = x1_pseudo / r_pseudo,
         w2_pseudo = x2_pseudo / r_pseudo)
# determine ideal value of k using the above
gw_pseudo <- gauss_gauge(pseudo_pred$w1_pseudo, pseudo_pred$w2_pseudo, test_radial$dep[2])
poss_k <- pseudo_pred$r_pseudo * gw_pseudo / ctau_fitted
k <- max(round(min(poss_k), 1) - 0.5, 1)

w_wts <- imp_weights(k, new_w, test_radial$dep[2], test_radial$alpha[2], "gauss")
r0w <- w_wts$r0w
rate <- w_wts$rate
w_wts <- w_wts$wts
resample_idx <- sample(1:length(new_w), replace = TRUE, size = 5000, prob = w_wts)
wstar <- new_w[resample_idx]
r0wstar <- k * r0w[resample_idx]
rstar <- qgamma(1 - runif(length(wstar)) * pgamma(r0wstar, shape = test_radial$alpha[2], 
                                                  rate = rate[resample_idx], lower.tail = F),
                shape = test_radial$alpha[2], rate = rate[resample_idx])
plot(x, y, pch = 20, xlim = c(0,15), ylim = c(0,15))
points(rstar * wstar, rstar * (1 - wstar), pch = 20, col = "red")
xstar <- cbind(rstar * wstar, rstar * (1 - wstar)) |> as_tibble() |> rename(x = V1, y = V2)

prob_x_given_r <- (which(
  with(
    xstar, 
    between(x, dim1[1], dim1[2]) & between(y, dim2[1], dim2[2]))
) |> length()) / nrow(xstar)

prob_r_over_thres <- length(data$idx)/data$N
prob_r_over_k_over_thres <- mean(w_wts)
prob_r_over_k <- prob_r_over_k_over_thres * prob_r_over_thres

prob_x_given_r * prob_r_over_k

rstar2 <- qgamma(1 - runif(length(new_w)) * pgamma(k * r0w, shape = test_radial$alpha[2], 
                                                  rate = rate, lower.tail = F),
                shape = test_radial$alpha[2], rate = rate)
points(rstar2 * new_w, rstar2 * (1 - new_w), pch = 20, col = "orange")
xstar2 <- cbind(rstar2 * new_w, rstar2 * (1 - new_w)) |> as_tibble() |> rename(x = V1, y = V2)

prob_x_given_r <- (which(
  with(
    xstar2, 
    between(x, dim1[1], dim1[2]) & between(y, dim2[1], dim2[2]))
) |> length()) / nrow(xstar2)



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
  cdf_eval <- pgamma(k * r0w, post_radial_alpha, gw)
  u_trunc <- runif(length(sim_w)) * (1-cdf_eval) + cdf_eval
  sim_r <- qgamma(u_trunc, post_radial_alpha, gw)
  qgamma(1 - runif(nsim) * pgamma(r0w_resample, shape = pars[1], 
                                  rate = rate_star, lower.tail = F),
         shape = pars[1], rate = rate_star)
  xy <- cbind(sim_r * sim_w, sim_r * (1 - sim_w)) |> as_tibble() |> rename(x = V1, y = V2)
  r_over_k <- pgamma(k * r0w, post_radial_alpha, gw, lower.tail = FALSE) / pgamma(r0w, post_radial_alpha, gw, lower.tail = FALSE)
  
  return(list(xy = xy, r0w = r0w, gw = gw, r_over_k = r_over_k))
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
  
  prob_pred <- length(in_box) / nrow(sim_data$xy) * mean(sim_data$r_over_k) * length(data_obs$idx)/data_obs$N
  
  return(list(prob_pred = prob_pred, in_box = in_box, k = k, xy = sim_data$xy))
}

preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, box, k_var, p) {
  post_radial <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_%s_radial.qs",
                               gauge, dep_type, dep_level, likelihood))
  post_ang_mix <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_ang_mix.qs",
                                dep_type, dep_level))
  post_ang_star <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_%s_ang_star.qs",
                                 dep_type, dep_level, gauge))
  
  mix <- map_dbl(1:200, function(i) {
    prob <- is_prob_pred(data_obs = fload(sprintf("data/%s/%s_%s.json", dep_type, dep_level, i)),
                         post_radial_params = post_radial[i, ],
                         post_angular_params = post_ang_mix[i, ],
                         box = box,
                         gauge_type = gauge, 
                         mix = TRUE,
                         k_var = k_var)$prob_pred
    p()
    prob
  })
  
  star <- map_dbl(1:200, function(i) {
    prob <- is_prob_pred(data_obs = fload(sprintf("data/%s/%s_%s.json", dep_type, dep_level, i)),
                         post_radial_params = post_radial[i, ],
                         post_angular_params = post_ang_star[i, ],
                         box = box,
                         gauge_type = gauge, 
                         mix = FALSE,
                         k_var = k_var)$prob_pred
    p()
    prob
  })
  
  preds <- cbind(mix, star) |> as_tibble() |>
    mutate(method = gauge,
           dataset = 1:200)
  return(preds)
}

preds_by_dep_level_lhood <- function(dep_type, dep_level, likelihood, box, k_var, p) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(lapply(gauge_library, function(x) preds_by_gauge(gauge = x, 
                                                          dep_type = dep_type, 
                                                          dep_level = dep_level, 
                                                          likelihood = likelihood, 
                                                          box = box,
                                                          k_var = k_var,
                                                          p = p)) |> 
           bind_rows())
}

weighted_preds_by_lhood <- function(dep_type, dep_level, likelihood, box, k_var, p) {
  preds <- preds_by_dep_level_lhood(dep_type = dep_type, 
                                    dep_level = dep_level, 
                                    likelihood = likelihood, 
                                    box = box, 
                                    k_var = k_var,
                                    p = p) |>
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

weighted_preds_by_level <- function(dep_type, dep_level, box, k_var, p) {
  lhood <- c("trunc", "cens")
  all_wts <- lapply(lhood, function(x) weighted_preds_by_lhood(dep_type = dep_type, 
                                                               dep_level = dep_level, 
                                                               likelihood = x,
                                                               box = box,
                                                               k_var = k_var,
                                                               p = p)) |> bind_rows()
  return(all_wts)
}

# create function that wraps everything to gether to make boxplot of predictions
create_predictions_boxplot <- function(dep_type, dep_level, box, k_var, p) {
  
  dim1 <- c(10, 12)
  if(box == "b1") {
    dim2 <- dim1
  } else if(box == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  plot_filename <- sprintf("figures/is_preds_boxplots/joint/%s_%s_%s_sim_angles.pdf",
                           dep_type, dep_level, box)
  qs_filename <-  sprintf("figures/is_preds_boxplots/joint/plot_objects/%s_%s_%s_sim_angles.qs",
                          dep_type, dep_level, box)
  plot_title <- sprintf("%s, %s, (%s) x (%s)",
                        dep_type, dep_level, paste(dim1, collapse = ","), paste(dim2, collapse = ","))
  
  preds_tib <- weighted_preds_by_level(dep_type = dep_type, 
                                       dep_level = dep_level, 
                                       box = box,
                                       k_var = k_var,
                                       p = p)
  
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
  p <- progressor(steps = nrow(all_combos) * 6 * 2 * 400)
  
  # Apply the function using apply and update the progress bar
  apply(all_combos, 1, function(row) {
    # p()  # Update the progress bar
    create_predictions_boxplot(dep_type = row["dep_types"],
                               dep_level = row["dep_levels"],
                               box = row["boxes"],
                               k_var = TRUE,
                               p = p)
  })
})
