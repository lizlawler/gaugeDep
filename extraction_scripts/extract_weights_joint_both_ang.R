# =============================================================================
# Computes BMA model weights across all 12 model variants (6 gauge functions
# x 2 angular density models: star and mix) for a single simulation study
# dataset. This is the "both angular" version of the weight extraction, which
# treats the angular density choice as part of the model comparison rather
# than fixing it. Run once per dataset via shell parallelism.
#
# Called by: shell_scripts/run_joint_model_wts.sh (one job per dataset)
# Inputs:    data/{dep_type}/{dep_level}_{data_num}.json
#            samplers/nimble/ang_mix_mcmc_fits/{dep_type}/{dep_level}_{data_num}.qs
#            samplers/rcpp/ang_star_mcmc_fits/{dep_type}/{gauge}_{dep_level}_{data_num}.qs
#            samplers/rcpp/radial_mcmc_fits/{dep_type}/{gauge}_{likelihood}_{dep_level}_{data_num}.qs
# Outputs:   fits_and_weights/wts_joint_model/both_ang/{dep_type}_{dep_level}_{likelihood}_{data_num}.qs
#
# Command-line args:
#   1. dep_type  -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level -- dependence strength ("low", "mid", "high")
#   3. data_num  -- dataset index (1–200)
#
# Note: currently only the censored likelihood is computed (the truncated
#   block is retained but commented out for reference).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type  <- args[1]
dep_level <- args[2]
data_num  <- args[3]

library(qs)
library(posterior)
library(loo)
library(dplyr)
library(tidyr)
library(gaugeDependence)

options(mc.cores = parallel::detectCores())

# Load data and pre-compute angular log-likelihoods --------------------------
datafile <- sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
data     <- RcppSimdJson::fload(datafile)
w        <- data$W

# Shared helpers for mixture density (avoids duplicating mix_lpdf / angular_loglik)
source("samplers/nimble/angular_mix_helpers.R")

# mix_ang_calc ---------------------------------------------------------------
# Computes the (n_iter x n_obs) pointwise loglik matrix for the angular
# mixture model on a single dataset.
mix_ang_calc <- function(angles, dep_type, dep_level, data_num) {
  params_file      <- sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
                              dep_type, dep_level, data_num)
  posterior_params <- qread(params_file)
  angular_loglik(angles, posterior_params)
}

mix_pw_loglik <- mix_ang_calc(angles = w, dep_type = dep_type,
                              dep_level = dep_level, data_num = data_num)

# star_ang_calc --------------------------------------------------------------
# Computes the (n_iter x n_obs) pointwise loglik matrix for the star-shaped
# angular density for a given gauge function.
star_ang_calc <- function(w, gauge, dep_type, dep_level, data_num) {
  params_file <- sprintf("samplers/rcpp/ang_star_mcmc_fits/%s/%s_%s_%s.qs",
                         dep_type, gauge, dep_level, data_num)
  params <- qread(params_file)$samples
  params <- params[, 1:(ncol(params) - 1)]  # drop sampler bookkeeping column
  angular_loglik(angles = w, dim = 2, posterior_params = params, gauge_type = gauge)
}

# Compute star angular logliks for all 6 gauges; store as named global objects
gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "rectangular", "dirichlet")
local({
  res <- setNames(
    lapply(gauge_library, function(x)
      star_ang_calc(w, x, dep_type, dep_level, data_num)),
    paste0("star_", gauge_library, "_pw_loglik")
  )
  list2env(res, envir = .GlobalEnv)
})
gc()

# radial_calc ----------------------------------------------------------------
# Computes the (n_iter x n_obs) pointwise loglik matrix for the radial model.
# Handles both truncated and censored likelihoods.
radial_calc <- function(data, gauge, likelihood, dep_type, dep_level, data_num) {
  params_file <- sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs",
                         dep_type, gauge, likelihood, dep_level, data_num)
  params <- qread(params_file)$samples
  params <- params[, 1:(ncol(params) - 1)]

  w   <- data$W
  r   <- data$R
  r0w <- data$r0_w
  idx <- data$idx

  if (likelihood == "trunc") {
    w <- w[idx]; r <- r[idx]; r0w <- r0w[idx]
  }

  radial_loglik(radii = r, threshold = r0w, angles = w,
                posterior_params = params,
                likelihood_type  = likelihood,
                gauge_type       = gauge)
}

# create_joint_loglik --------------------------------------------------------
# Adds the radial and angular log-likelihood matrices for a given gauge and
# likelihood type. For "star" angle_dens, uses the gauge-matched star loglik;
# for "mix", uses the single gauge-agnostic mixture loglik.
create_joint_loglik <- function(data, gauge, likelihood, angle_dens) {
  temp_radial  <- get(sprintf("radial_%s_%s_pw_loglik", gauge, likelihood))
  temp_angular <- if (angle_dens == "star") {
    get(sprintf("star_%s_pw_loglik", gauge))
  } else {
    get("mix_pw_loglik")
  }

  if (likelihood == "trunc") {
    temp_angular <- temp_angular[, data$idx]
  }
  temp_radial + temp_angular
}

# extract_lpd_pt / create_lpd_list / model_weights ---------------------------
# Standard LOO-ELPD extraction and weight computation across all 12 models.

extract_lpd_pt <- function(data, gauge, likelihood, angle_dens) {
  temp     <- create_joint_loglik(data, gauge, likelihood, angle_dens)
  loo_temp <- loo(temp)
  loo_temp$pointwise[, "elpd_loo"]
}

create_lpd_list <- function(data, likelihood) {
  star_lpd <- sapply(gauge_library,
                     function(x) extract_lpd_pt(data, x, likelihood, "star"))
  colnames(star_lpd) <- paste0(gauge_library, "_star")

  mix_lpd <- sapply(gauge_library,
                    function(x) extract_lpd_pt(data, x, likelihood, "mix"))
  colnames(mix_lpd) <- paste0(gauge_library, "_mix")

  cbind(star_lpd, mix_lpd)
}

model_weights <- function(data, likelihood) {
  temp <- create_lpd_list(data, likelihood)
  list(
    stacking         = stacking_weights(temp),
    pseudobma_boot   = pseudobma_weights(temp),
    pseudobma_noboot = pseudobma_weights(temp, BB = FALSE)
  )
}

# Compute and save censored-likelihood weights --------------------------------
# (Truncated block retained below for reference but currently inactive)
local({
  res <- setNames(
    lapply(gauge_library, function(x)
      radial_calc(data, x, "cens", dep_type, dep_level, data_num)),
    paste0("radial_", gauge_library, "_cens_pw_loglik")
  )
  list2env(res, envir = .GlobalEnv)
})
gc()

qsave(model_weights(data, "cens"),
      sprintf("fits_and_weights/wts_joint_model/both_ang/%s_%s_%s_%s.qs",
              dep_type, dep_level, "cens", data_num))

# Truncated likelihood (inactive — uncomment to enable):
# local({
#   res <- setNames(
#     lapply(gauge_library, function(x)
#       radial_calc(data, x, "trunc", dep_type, dep_level, data_num)),
#     paste0("radial_", gauge_library, "_trunc_pw_loglik")
#   )
#   list2env(res, envir = .GlobalEnv)
# })
# gc()
# qsave(model_weights(data, "trunc"),
#       sprintf("fits_and_weights/wts_joint_model/both_ang/%s_%s_%s_%s.qs",
#               dep_type, dep_level, "trunc", data_num))
