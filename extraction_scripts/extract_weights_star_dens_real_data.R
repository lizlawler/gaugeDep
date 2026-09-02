# =============================================================================
# Computes BMA model weights using only the star-shaped angular density
# log-likelihoods for both real data stations (Friend Mountain and Redstone).
# This is the angular-component-only weight extraction, used to compare gauge
# functions on the basis of how well they model the angular distribution alone.
#
# Run interactively (no command-line args).
# Inputs:    samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/{data_type}_{gauge}.qs
# Outputs:   fits_and_weights/wts_joint_model/{data_type}_star_dens_only.qs
# =============================================================================

library(qs)
library(posterior)
library(loo)
library(dplyr)
library(tidyr)

options(mc.cores = parallel::detectCores())

# extract_lpd_pt -------------------------------------------------------------
# Runs LOO-CV on the angular pointwise loglik for one station and gauge.
extract_lpd_pt <- function(data_type, gauge) {
  temp     <- qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/pw_loglik/%s_%s.qs",
                            data_type, gauge))
  loo_temp <- loo(temp)
  loo_temp$pointwise[, "elpd_loo"]
}

# create_lpd_list ------------------------------------------------------------
# Returns a named list of ELPD vectors across all 6 gauges for one station.
create_lpd_list <- function(data_type) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  setNames(
    sapply(gauge_library, function(x) extract_lpd_pt(x, data_type = data_type)),
    gauge_library
  )
}

# model_weights --------------------------------------------------------------
# Computes and saves all three BMA weight variants for one station.
model_weights <- function(data_type) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  temp <- create_lpd_list(data_type)

  wts <- tibble(
    stacking         = stacking_weights(temp),
    pseudobma_boot   = pseudobma_weights(temp),
    pseudobma_noboot = pseudobma_weights(temp, BB = FALSE)
  ) |>
    mutate(method            = gauge_library,
           stacking          = as.numeric(stacking),
           pseudobma_boot    = as.numeric(pseudobma_boot),
           pseudobma_noboot  = as.numeric(pseudobma_noboot))

  wts_file <- sprintf("fits_and_weights/wts_joint_model/%s_star_dens_only.qs", data_type)
  qsave(wts, wts_file)
}

# Run for both stations
model_weights("redstone")
model_weights("friendmtn")
