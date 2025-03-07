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
  mc_star <- est_volume(n = 100, pars, gauge_type)
  gauge_fn <- get_gauge_function(gauge_type)
  gw <- gauge_fn(w1, w2, pars)
  return(1 / (gw^2 * 2 * mc_star))
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

true_inv_log_prob <- function(dim1, dim2, dep) {
  dim1_star <- 1 / dim1
  dim2_star <- 1 / dim2
  
  margins <- c(1, 1, 1)
  upper_right <- pbvevd(q = c(dim1_star[2], dim2_star[2]), model = "log", dep = dep, mar1 = margins)
  upper_left <- pbvevd(q = c(dim1_star[1], dim2_star[2]), model = "log", dep = dep, mar1 = margins)
  lower_right <- pbvevd(q = c(dim1_star[2], dim2_star[1]), model = "log", dep = dep, mar1 = margins)
  lower_left <- pbvevd(q = c(dim1_star[1], dim2_star[1]), model = "log", dep = dep, mar1 = margins)
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
    filter(x1 >= 0, x2 >= 0) |> # for B3, some samples may be < 0; but the probability is next to zero, so don't need to account for in density eval later on
    mutate(r = x1 + x2, w1 = x1 / r, w2 = x2 / r)
  
  return(is_samp_mvn)
}

# probability calculation
is_prob_pred <- function(imp_samples = NULL, 
                         post_radial_params, 
                         post_angular_params,
                         box = "b1",
                         gauge_type, 
                         ang_dens = "star",
                         sir = FALSE) {
  
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
  
  # estimate RW density
  gauge_vals <- gauge_fn(imp_samples$w1, imp_samples$w2, post_radial_dep)
  r0w <- qgamma(0.95, shape = post_radial_alpha, rate = gauge_vals)
  r_giv_w_dens <- trunc_gamma(imp_samples$r, r0w, alpha = post_radial_alpha, beta = gauge_vals)
  
  rw_dens <- r_giv_w_dens * ang_dens
  is_dens <- mvtnorm::dmvnorm(imp_samples[, 1:2], mean = c(mean(dim1), mean(dim2)), sigma = 2 * diag(2)) * imp_samples$r
  wts <- rw_dens/is_dens
  
  # resample points using importance weights (sampling importance resampling)
  if(sir) {
    sir_idx <- sample(1:nrow(imp_samples), size = floor(nrow(imp_samples)/15), 
                      replace = FALSE, prob = wts)
    imp_samples <- imp_samples[sir_idx,]
    wts <- wts[sir_idx,]
  }
  
  # indicator of being in box
  idx_in_box <- which(
    with(
      imp_samples,
      between(r, dim1[1] / w1, dim1[2] / w1) & 
        between(r, dim2[1] / (1-w1), dim2[2] / (1-w1)))
  )
  
  return(sum(wts[idx_in_box]) / nrow(imp_samples) * 0.05)
}
# 
# test_radial <- qread("fits_and_weights/post_params_joint/gauss_gauss_high_wc_trunc_radial.qs")
# test_ang_mix <- qread("fits_and_weights/post_params_joint/gauss_high_wc_ang_mix.qs")
# test_ang_star <- qread("fits_and_weights/post_params_joint/gauss_high_wc_gauss_ang_star.qs")
# 
# box <- "b3"
# dim1 <- c(10, 12)
# dim2 <- case_when(
#   box == "b1" ~ dim1,
#   box == "b2" ~ c(6, 8),
#   TRUE ~ c(2, 4)
# )

# is_samps <- gen_is_samples(box,total_n = 50000)
# data_num <- 40
# data <- RcppSimdJson::fload(sprintf("data/gauss/high_wc_%s.json", data_num))
# r <- data$R
# w <- data$W
# ang_dens <- mix_dens(is_samps$w1, test_ang_mix[data_num,])
# gw <- gauss_gauge(is_samps$w1, is_samps$w2, test_radial$dep[data_num])
# r0w <- qgamma(0.95, test_radial$alpha[data_num], gw)
# r_giv_w <- trunc_gamma(is_samps$r, r0w, test_radial$alpha[data_num], gw)
# rw_dens <- r_giv_w * ang_dens
# is_dens <- mvtnorm::dmvnorm(is_samps[, 1:2], mean = c(mean(dim1), mean(dim2)), sigma = 2 * diag(2)) * is_samps$r
# wts <- rw_dens / is_dens
# 
# sir_idx <- sample(1:nrow(is_samps),size = 3333, replace = FALSE, prob = wts)
# plot(r * w, r * (1-w), pch = 20, xlim = c(0,16), ylim = c(0,16))
# points(is_samps$x1, is_samps$x2,pch = 20, col = "blue")
# points(is_samps$x1[sir_idx], is_samps$x2[sir_idx], col = "red", pch = 20)
# idx_in_box<- which(
#   with(
#     is_samps[sir_idx,],
#     between(r, dim1[1] / w1, dim1[2] / w1) &
#       between(r, dim2[1] / (1-w1), dim2[2] / (1-w1)))
# )
# 
# points(is_samps$x1[sir_idx][idx_in_box], is_samps$x2[sir_idx][idx_in_box], col = "orange", pch = 20)
# # 
# sum(wts[sir_idx][idx_in_box]) / length(sir_idx) * 0.05
# true_gauss_prob(dim1, dim2, 0.8)
# 
# idx_in_box_is <- which(
#   with(
#     is_samps,
#     between(r, dim1[1] / w1, dim1[2] / w1) &
#       between(r, dim2[1] / (1-w1), dim2[2] / (1-w1)))
# )
# 
# sum(wts[idx_in_box_is]) / length(wts) * 0.05

preds_by_gauge <- function(gauge, dep_type, dep_level, likelihood, box, sir = FALSE) {
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
                 ang_dens = "mix",
                 sir = sir)
  })
  
  star <- map_dbl(1:200, function(i) {
    is_prob_pred(post_radial_params = post_radial[i, ],
                 post_angular_params = post_ang_star[i, ],
                 box = box,
                 gauge_type = gauge,
                 ang_dens = "star",
                 sir = sir)
  })
  
  preds <- cbind(mix, star) |> as_tibble() |>
    mutate(method = gauge,
           dataset = 1:200)
  return(preds)
}

# test <- preds_by_gauge("gauss", "gauss", "high", "trunc", "b1")

preds_by_dep_level_lhood <- function(dep_type, dep_level, likelihood, box, sir = FALSE) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(lapply(gauge_library, function(x) preds_by_gauge(gauge = x, 
                                                          dep_type = dep_type, 
                                                          dep_level = dep_level, 
                                                          likelihood = likelihood, 
                                                          box = box,
                                                          sir = sir)) |> 
           bind_rows())
}

weighted_preds_by_lhood <- function(dep_type, dep_level, likelihood, box, sir = FALSE) {

  preds <- preds_by_dep_level_lhood(dep_type = dep_type, 
                                    dep_level = dep_level, 
                                    likelihood = likelihood, 
                                    box = box, 
                                    sir = sir) |>
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

weighted_preds_by_level <- function(dep_type, dep_level, box, sir = FALSE) {
  lhood <- c("trunc", "cens")
  all_wts <- lapply(lhood, function(x) weighted_preds_by_lhood(dep_type = dep_type, 
                                                               dep_level = dep_level, 
                                                               likelihood = x,
                                                               box = box,
                                                               sir = sir)) |> bind_rows()
  return(all_wts)
}

# create function that wraps everything to gether to make boxplot of predictions
create_predictions_boxplot <- function(dep_type, dep_level, box, sir = FALSE) {
  
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
                                       box = box,
                                       sir = sir)
  
  # determine true probability
  if(dep_type == "gauss") {
    levels_list <- list(low = 0.1, mid = 0.5, high = 0.9, high_wc = 0.8)
    true_prob <- true_gauss_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]))
  } else if(dep_type == "logistic"){
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1, low_wc = 0.8, mid_wc = 0.4)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "log")
  } else {
    levels_list <- list(low = 0.1, mid = 1, high = 3)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "hr")
  }
  
  # create boxplot
  plot <- preds_tib |> ggplot(aes(x = ang_dens, y = preds, fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = true_prob, col = "darkgrey", linetype = "longdash") +
    ggtitle(plot_title) + 
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
  ggsave(plot_filename,
         plot = plot,
         bg = 'transparent',
         width = 8,
         height = 7,
         dpi = 320)
  qsave(plot, qs_filename)
  print(paste0(plot_filename, " has been saved"))
}

create_predictions_boxplot("gauss","mid", "b1")
create_predictions_boxplot("gauss","mid", "b2")
create_predictions_boxplot("gauss","mid", "b3")

dep_types <- c("gauss", "logistic")
dep_levels <- c("high_wc", "mid_wc", "low_wc")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)
all_combos <- all_combos |> filter(!(dep_types == "gauss" & dep_levels != "high_wc"),
                                   !(dep_types == "logistic" & dep_levels == "high_wc"))

dep_types <- c("husler_reiss")
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
                               box = row["boxes"])
    p()  # Update the progress bar
  })
})
