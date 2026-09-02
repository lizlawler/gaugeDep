# =============================================================================
# Computes BMA-weighted importance-sampling predictions using the "both angular"
# weight scheme, in which stacking weights are computed jointly across all 12
# model variants (6 gauge functions x 2 angular density models). Run after
# merge_wts_both_ang.R has produced the combined weight files.
#
# Inputs:    fits_and_weights/post_params_joint/...qs
#            fits_and_weights/wts_joint_model/{dep_type}_{dep_level}_{likelihood}_both_ang.qs
# Outputs:   figures/is_preds_boxplots/joint/pred_tibbles/both_ang/{dep_type}_{dep_level}_{likelihood}_{box}.qs
#            figures/is_preds_boxplots/joint/mse_tables/both_ang/{...}.qs
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

make_wts_df <- function(dep_type, dep_level, likelihood) {
  comb_wts_file <- sprintf("fits_and_weights/wts_joint_model/%s_%s_%s_both_ang.qs",
                      dep_type, dep_level, likelihood)
  
  # Check if BMA weights file already exists
  if (file.exists(comb_wts_file)) {
    return(qread(comb_wts_file))
  }
  
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  ang_lib <- c("star", "mix")
  gauge_ang_lib <- paste(gauge_library, rep(ang_lib, times = c(6,6)), sep = "_")
  
  wts_files <- list.files(path = "fits_and_weights/wts_joint_model/both_ang/",
                          pattern = sprintf("%s_%s_%s", dep_type, dep_level, likelihood),
                          full.names = TRUE)
  combined_wts <- map_dfr(wts_files, function(f) {
    df <- qread(f) |> bind_cols() |> mutate(method = gauge_ang_lib)
    df |> mutate(dataset = as.integer(str_extract(basename(f), "\\d+"))) |>
      mutate(stacking = as.numeric(stacking),
             pseudobma_boot = as.numeric(pseudobma_boot),
             pseudobma_noboot = as.numeric(pseudobma_noboot))
  })
  qsave(combined_wts, comb_wts_file)
  combined_wts
}

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
    filter(x1 >= 0, x2 >= 0) |> # adjust for this "truncation" in the is_prob_pred function
    mutate(r = x1 + x2, w1 = x1 / r, w2 = x2 / r)
  
  return(is_samp_mvn)
}

# probability calculation
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
  
  # estimate RW density
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
           dataset = as.integer(1:200))
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
    pivot_longer(cols = c(mix, star), names_to = "ang_dens", values_to = "preds") |>
    mutate(method = paste0(method, "_", ang_dens)) |> select(-ang_dens)
  wts <- make_wts_df(dep_type = dep_type, dep_level = dep_level, likelihood = likelihood)
  wtd_preds <- suppressMessages(preds |> left_join(wts) |>
                                  mutate(stacking_preds = preds * stacking,
                                         pseudo_boot = pseudobma_boot * preds,
                                         pseudo_noboot = pseudobma_noboot * preds) |>
                                  group_by(dataset) |>
                                  summarize(stacking_predictions = sum(stacking_preds),
                                            pseudobma_boot_preds = sum(pseudo_boot),
                                            pseudobma_noboot_preds = sum(pseudo_noboot)) |>
                                  ungroup())
  boxplot_wts <- wtd_preds |>
    pivot_longer(cols = -dataset, names_to = "method", values_to = "preds") |>
    mutate(method = case_when(grepl("stacking", method) ~ 'Stacking',
                              grepl("noboot", method) ~ 'Pseudo-BMA',
                              grepl("boot", method) ~ 'Pseudo-BMA+'),
           method = as.factor(method),
           model = likelihood)
  return(boxplot_wts)
}

# weighted_preds_by_level <- function(dep_type, dep_level, box) {
#   lhood <- c("trunc", "cens")
#   all_wts <- lapply(lhood, function(x) weighted_preds_by_lhood(dep_type = dep_type, 
#                                                                dep_level = dep_level, 
#                                                                likelihood = x,
#                                                                box = box)) |> bind_rows()
#   return(all_wts)
# }

# create function that wraps everything to gether to make boxplot of predictions
create_predictions_boxplot <- function(dep_type, dep_level, likelihood, box) {
  
  dim1 <- c(10, 12)
  if(box == "b1") {
    dim2 <- dim1
  } else if(box == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  preds_tib <- weighted_preds_by_lhood(dep_type = dep_type,
                                       dep_level = dep_level,
                                       likelihood = likelihood,
                                       box = box)
  qsave(preds_tib, sprintf("figures/is_preds_boxplots/joint/pred_tibbles/both_ang/%s_%s_%s_%s.qs",
                           dep_type, dep_level, likelihood, box))
  
  true_prob <- get_true_prob(dep_type, dep_level, dim1, dim2)
  
  # create boxplot
  plot <- preds_tib |> ggplot(aes(x = model, y = pmax(preds, .Machine$double.eps), fill = method)) +
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
    # scale_x_discrete(labels=c("Cens., Mix", "Cens., Star", "Trunc., Mix", "Trunc., Star")) +
    xlab("Likelihood") + ylab("Prediction probabilities")
  ggsave(sprintf("figures/is_preds_boxplots/joint/both_ang/%s_%s_%s_%s.pdf", dep_type, dep_level, likelihood, box),
         plot = plot,
         bg = 'transparent',
         width = 8,
         height = 7,
         dpi = 320)
  qsave(plot, sprintf("figures/is_preds_boxplots/joint/plot_objects/both_ang/%s_%s_%s_%s.qs",dep_type, dep_level, likelihood, box))
  print(sprintf("Boxplot for %s %s %s, in box %s has been saved", dep_level, dep_type, likelihood, box))

  mse_table <- preds_tib |> 
    mutate(truth = true_prob,
           diff = preds - truth) |>
    group_by(method) |>
    summarise(mse = mean(diff^2)) |>
    ungroup() |>
    mutate(log_mse = log(mse), rmse_norm = sqrt(mse)/true_prob) |>
    arrange(rmse_norm)
  qsave(mse_table, sprintf("figures/is_preds_boxplots/joint/mse_tables/both_ang/%s_%s_%s_%s.qs",dep_type, dep_level, likelihood, box))
}

dep_types <- c("gauss", "logistic")
dep_levels <- c("high", "mid", "low")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)

with_progress({
  p <- progressor(steps = nrow(all_combos))
  
  future_pmap(all_combos, function(dep_types, dep_levels, boxes) {
    create_predictions_boxplot(dep_type = dep_types,
                               dep_level = dep_levels,
                               likelihood = "cens",
                               box = boxes)
    p()  # Update the progress bar
  }, .options = furrr_options(seed = TRUE))
})

