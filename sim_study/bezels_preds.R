# =============================================================================
# Computes bivariate exceedance probability predictions using the BezELS
# (Bezier extreme level set) radial model combined with the angular mixture
# (NIMBLE stick-breaking) density, via importance sampling. BezELS provides
# a nonparametric alternative to the gauge-function radial model and serves
# as a competitor method in the simulation study.
#
# Run after bezels_mcmc.R and the angular mixture parameter extraction have
# completed. Parallelises across scenarios via furrr.
#
# Inputs:    samplers/bezels/radial_bezels_fits/{dep_type}/{dep_level}_{i}.qs
#            fits_and_weights/post_params_joint/{dep_type}_{dep_level}_ang_mix.qs
# Outputs:   figures/is_preds_boxplots/joint/pred_tibbles/{dep_type}_{dep_level}_{box}_bezels.qs
#            figures/is_preds_boxplots/joint/mse_tables/{dep_type}_{dep_level}_{box}_bezels.qs
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(progressr) 
library(RcppSimdJson)
library(gaugeDependence)
library(qs)
library(evd)
library(grafify)
library(stringr)
library(BezELS)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

# Shared ground-truth probability functions
source("sim_study/true_prob_utils.R")

library(furrr)
plan(multisession, workers = parallel::detectCores()/2)

trunc_gamma <- function(x, xmin, alpha, beta) {
  unnorm_pdf <- dgamma(x, shape = alpha, rate = beta, log = TRUE)
  norm_cst <- pgamma(xmin, shape = alpha, rate = beta, lower.tail = F, log.p = TRUE)
  return(exp(unnorm_pdf - norm_cst))
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

post_gamma_rate <- function(mcmc_samples, w_mat, thin.by = 10) {
  n.data = nrow(w_mat)
  n.samples = nrow(mcmc_samples)
  sss = seq(1, n.samples, by = thin.by)
  
  n.post = length(sss)
  post_rate <- matrix(NA, nrow = n.post, ncol = n.data)
  
  p0y <- (mcmc_samples[sss, 1])
  p1x <- (mcmc_samples[sss, 2])
  p1y <- (mcmc_samples[sss, 3])
  p2x <- (mcmc_samples[sss, 4])
  p3 <- (mcmc_samples[sss, 5])
  p4y <- (mcmc_samples[sss, 6])
  p5x <- (mcmc_samples[sss, 7])
  p5y <- (mcmc_samples[sss, 8])
  p6x <- (mcmc_samples[sss, 9])
  
  for (i in 1:n.post) {
    p <- matrix(c(0, p0y[i], p1x[i], p1y[i], p2x[i], 1, p3[i], 
                  p3[i], 1, p4y[i], p5x[i], p5y[i], p6x[i], 0), 
                nrow = 7, ncol = 2, byrow = TRUE)
    post_rate[i, ] <- gx(N = n.data, p, w_mat[, 2]/w_mat[, 1], w_mat[, 1])
  }
  
  return(apply(post_rate, 2, median))
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
                         post_bezels_alpha, 
                         bezels_samples,
                         post_angular_params,
                         box = "b1") {
  
  if(is.null(imp_samples)) {
    imp_samples <- gen_is_samples(box = box)
  }
  
  # specify dimensions of box
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  # calculate posterior gamma rate using bezels samples
  angles <- cbind(imp_samples$w1, imp_samples$w2)
  post_bezels_rate <- post_gamma_rate(mcmc_samples = bezels_samples, w_mat = angles, thin.by = 20)
  
  # compute angular density
  ang_dens <- mix_dens(imp_samples$w1, post_angular_params)
  
  # estimate RW density
  r0w <- qgamma(0.95, shape = post_bezels_alpha, rate = post_bezels_rate)
  r_giv_w_dens <- trunc_gamma(imp_samples$r, r0w, alpha = post_bezels_alpha, beta = post_bezels_rate)
  
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

preds_by_dataset <- function(dep_type, dep_level, dataset, box) {
  bezels_mcmc <- qread(sprintf("samplers/bezels/radial_bezels_fits/%s/%s_%s.qs", dep_type, dep_level, dataset))
  post_alpha <- median(bezels_mcmc[,10])
  post_ang_mix <- qread(sprintf("fits_and_weights/post_params_joint/%s_%s_ang_mix.qs",
                                dep_type, dep_level))[dataset,]
  
  return(is_prob_pred(post_bezels_alpha = post_alpha,
                      bezels_samples = bezels_mcmc,
                      post_angular_params = post_ang_mix,
                      box = box))
}

preds_by_box <- function(dep_type, dep_level, box, p) {
  
  preds <- map_dbl(1:200, function(i) {
    prob <- preds_by_dataset(dep_type = dep_type,
                             dep_level = dep_level,
                             dataset = i,
                             box = box)
    p()  # Update the progress bar
    prob
  })
  
  return(preds |> as_tibble() |> rename(preds = value) |> mutate(dataset = 1:200))
}

create_predictions_tibble <- function(dep_type, dep_level, box, p) {
  
  dim1 <- c(10, 12)
  if(box == "b1") {
    dim2 <- dim1
  } else if(box == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  preds_tib <- preds_by_box(dep_type = dep_type,
                            dep_level = dep_level,
                            box = box,
                            p = p) |>
    mutate(ang_dens = "BezELS", method = NA)
  qsave(preds_tib, sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_bezels.qs", dep_type, dep_level, box))
  
  true_prob <- get_true_prob(dep_type, dep_level, dim1, dim2)
  
  preds_tib <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_bezels.qs", dep_type, dep_level, box))
  
  mse_table <- preds_tib |> 
    mutate(truth = true_prob,
           diff = preds - truth) |>
    group_by(method, ang_dens) |>
    summarise(mse = mean(diff^2)) |>
    ungroup() |>
    mutate(log_mse = log(mse), rmse_norm = sqrt(mse)/true_prob)
  qsave(mse_table, sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_bezels.qs",dep_type, dep_level, box))
}

dep_types <- c("gauss", "logistic", "husler_reiss")
dep_levels <- c("high", "mid", "low")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)

with_progress({
  p <- progressor(steps = nrow(all_combos) * 200)
  
  future_pmap(all_combos, function(dep_types, dep_levels, boxes) {
    create_predictions_tibble(dep_type = dep_types,
                              dep_level = dep_levels,
                              box = boxes,
                              p = p)
  }, .options = furrr_options(seed = TRUE))
})

