# =============================================================================
# Extracts posterior parameters from real fire weather data MCMC fits for both
# stations (Friend Mountain, Redstone). Provides three extraction functions
# covering each model component: radial (Rcpp), angular star (Rcpp), and
# angular mixture (NIMBLE). Intended to be run interactively rather than
# via command line — set data_name at the top to select the station.
#
# Inputs:
#   samplers/rcpp/radial_mcmc_fits/real_data/{data}_{gauge}_{likelihood}_2chains.qs
#   samplers/rcpp/ang_star_mcmc_fits/real_data/{data}_{gauge}_2chains.qs
#   samplers/nimble/ang_mix_mcmc_fits/real_data/{data}_2chains.qs
# =============================================================================

library(qs)
library(dplyr)
library(tidyr)

data_name <- "friendmtn"  # set to "redstone" for the other station

# extract_post_params_radial -------------------------------------------------
# Returns posterior samples or summary (median) for the radial component.
# If summarize = TRUE, returns a one-row tibble of posterior medians.
# If summarize = FALSE, returns the full draws tibble.
extract_post_params_radial <- function(gauge, likelihood, data, summarize = TRUE) {
  params <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/real_data/%s_%s_%s_2chains.qs",
                          data, gauge, likelihood)) |>
    as_tibble()

  if (summarize) {
    params <- params |>
      apply(MARGIN = 2, FUN = median) |>
      t() |>
      as_tibble()
  }
  return(params)
}

# extract_post_params_ang_star -----------------------------------------------
# Returns posterior samples or summary for the star-shaped angular density.
extract_post_params_ang_star <- function(gauge, data, summarize = TRUE) {
  params <- qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/%s_%s_2chains.qs",
                          data, gauge)) |>
    as_tibble()

  if (summarize) {
    params <- params |>
      apply(MARGIN = 2, FUN = median) |>
      t() |>
      as_tibble()
  }
  return(params)
}

# extract_post_params_ang_mix ------------------------------------------------
# Returns posterior samples or per-chain means for the angular mixture model.
# If summarize = TRUE, returns per-chain posterior means (one row per chain)
# via group_by(.chain) to preserve chain-level information.
extract_post_params_ang_mix <- function(data, summarize = TRUE) {
  params <- qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/%s_2chains.qs",
                          data)) |>
    as_tibble() |>
    select(matches("probs|alphastar|betastar|chain_id"))

  if (summarize) {
    params <- params |>
      group_by(chain_id) |>
      summarize(across(where(is.numeric), mean), .groups = "drop")
  }
  return(params)
}
