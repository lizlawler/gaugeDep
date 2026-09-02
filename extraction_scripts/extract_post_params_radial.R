# =============================================================================
# Extracts and summarises posterior parameters from the Rcpp radial MCMC fits
# across all 200 simulation study datasets. Saves posterior medians of the
# radial dependence and shape parameters for downstream visualisation.
#
# Called by: shell_scripts/call_radial_params.sh
# Inputs:    samplers/rcpp/radial_mcmc_fits/{dep_type}/{gauge}_{likelihood}_{dep_level}_{i}.qs
# Outputs:   fits_and_weights/post_params_joint/{gauge}_{dep_type}_{dep_level}_{likelihood}_radial.qs
#
# Command-line args:
#   1. dep_type   -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level  -- dependence strength ("low", "mid", "high")
#   3. gauge      -- gauge function ("gauss", "logistic", "inv_log", etc.)
#   4. likelihood -- likelihood type ("trunc" or "cens")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type   <- args[1]
dep_level  <- args[2]
gauge      <- args[3]
likelihood <- args[4]

library(qs)
library(dplyr)
library(tidyr)

# extract_post_params_radial -------------------------------------------------
# Loads Rcpp MCMC output for one dataset and returns a one-row tibble of
# posterior medians for all radial parameters (excluding sigma_m).
extract_post_params_radial <- function(dep_type, dep_level, gauge, likelihood, data_num) {
  params <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs",
                          dep_type, gauge, likelihood, dep_level, data_num))$samples |>
    as_tibble()

  params |>
    select(-sigma_m) |>  # drop sampler bookkeeping column
    apply(MARGIN = 2, FUN = median) |>
    t() |>
    as_tibble() |>
    mutate(dataset = data_num)
}

# create_tib_post_params_radial ----------------------------------------------
# Loops over all 200 datasets, stacks the one-row summaries, and saves.
create_tib_post_params_radial <- function(dep_type, dep_level, gauge, likelihood) {
  tib_post_params <- sapply(
    1:200,
    function(x) extract_post_params_radial(dep_type, dep_level, gauge, likelihood, x),
    simplify = FALSE
  ) |>
    bind_rows()

  filepath <- sprintf("fits_and_weights/post_params_joint/%s_%s_%s_%s_radial.qs",
                      gauge, dep_type, dep_level, likelihood)
  qsave(tib_post_params, filepath)
  print("Posterior medians of radial parameters have been successfully saved")
}

create_tib_post_params_radial(dep_type = dep_type, dep_level = dep_level,
                              gauge = gauge, likelihood = likelihood)
