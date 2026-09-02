# =============================================================================
# Computes BMA model weights (stacking, pseudo-BMA, pseudo-BMA+) for the joint
# model on real fire weather data. Mirrors extract_weights_joint.R but for
# a single station rather than looping over simulated datasets.
#
# Called by: shell_scripts/local_machine/run_joint.sh (or similar)
# Inputs:    samplers/rcpp/radial_mcmc_fits/real_data/pw_loglik/...qs
#            samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/...qs (star)
#            samplers/nimble/ang_mix_mcmc_fits/real_data/pw_loglik/...qs (mix)
#            data/raw/{data_type}_expo.qs                                (for idx)
# Outputs:   fits_and_weights/wts_joint_model/{data_type}_{likelihood}_{angle_dens}.qs
#            samplers/joint_loglik/real_data/...qs                       (cached)
#
# Command-line args:
#   1. data_type  -- station identifier ("friendmtn" or "redstone")
#   2. likelihood -- likelihood type ("trunc" or "cens")
#   3. angle_dens -- angular density model ("star" or "mix")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
data_type  <- args[1]
likelihood <- args[2]
angle_dens <- args[3]

library(qs)
library(posterior)
library(loo)
library(dplyr)
library(tidyr)

options(mc.cores = parallel::detectCores())

# create_joint_loglik --------------------------------------------------------
# Assembles the joint pointwise log-likelihood matrix for one gauge function
# by summing the radial and angular component matrices. Caches to disk.
create_joint_loglik <- function(gauge, likelihood = "cens", angle_dens = "star") {
  loglik_file <- sprintf("samplers/joint_loglik/real_data/%s_%s_%s_%s.qs",
                         data_type, angle_dens, gauge, likelihood)

  if (file.exists(loglik_file)) return(qread(loglik_file))

  temp_radial <- qread(sprintf(
    "samplers/rcpp/radial_mcmc_fits/real_data/pw_loglik/%s_%s_%s.qs",
    data_type, gauge, likelihood
  ))

  temp_angular <- if (angle_dens == "star") {
    qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/%s_%s.qs",
                  data_type, gauge))
  } else {
    qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/pw_loglik/%s.qs",
                  data_type))
  }

  # For truncated likelihood, restrict angular loglik to exceedance indices only
  if (likelihood == "trunc") {
    idx          <- qread(sprintf("data/raw/%s_expo.qs", data_type))$idx
    temp_angular <- temp_angular[, idx]
  }

  temp_joint <- temp_radial + temp_angular
  qsave(temp_joint, loglik_file)
  temp_joint
}

# extract_lpd_pt / create_lpd_list / model_weights ---------------------------
# LOO-ELPD extraction and BMA weight computation across all 6 gauge functions.

extract_lpd_pt <- function(gauge, likelihood, angle_dens) {
  temp     <- create_joint_loglik(gauge, likelihood, angle_dens)
  loo_temp <- loo(temp)
  loo_temp$pointwise[, "elpd_loo"]
}

create_lpd_list <- function(likelihood, angle_dens) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  setNames(
    sapply(gauge_library,
           function(x) extract_lpd_pt(x, likelihood, angle_dens)),
    gauge_library
  )
}

model_weights <- function(likelihood, angle_dens) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  temp <- create_lpd_list(likelihood, angle_dens)

  wts <- tibble(
    stacking         = stacking_weights(temp),
    pseudobma_boot   = pseudobma_weights(temp),
    pseudobma_noboot = pseudobma_weights(temp, BB = FALSE)
  ) |>
    mutate(method            = gauge_library,
           stacking          = as.numeric(stacking),
           pseudobma_boot    = as.numeric(pseudobma_boot),
           pseudobma_noboot  = as.numeric(pseudobma_noboot))

  wts_file <- sprintf("fits_and_weights/wts_joint_model/%s_%s_%s.qs",
                      data_type, likelihood, angle_dens)
  qsave(wts, wts_file)
  print(sprintf("Weights extracted and saved for %s, angular density: %s, likelihood: %s",
                data_type, angle_dens, likelihood))
}

model_weights(likelihood = likelihood, angle_dens = angle_dens)
