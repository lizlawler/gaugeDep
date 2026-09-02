# =============================================================================
# Extracts and summarises posterior parameters from the Rcpp star-shaped
# angular density MCMC fits across all 200 simulation study datasets. Saves
# posterior medians of the angular dependence parameter for downstream
# visualisation and prediction.
#
# Called by: shell_scripts/call_ang_star_params.sh
# Inputs:    samplers/rcpp/ang_star_mcmc_fits/{dep_type}/{gauge}_{dep_level}_{i}.qs
# Outputs:   fits_and_weights/post_params_joint/{dep_type}_{dep_level}_{gauge}_ang_star.qs
#
# Command-line args:
#   1. dep_type  -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level -- dependence strength ("low", "mid", "high")
#   3. gauge     -- gauge function ("gauss", "logistic", "inv_log", etc.)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type  <- args[1]
dep_level <- args[2]
gauge     <- args[3]

library(qs)
library(dplyr)
library(tidyr)

# extract_post_params_ang_star -----------------------------------------------
# Loads Rcpp MCMC output for one dataset and returns a one-row tibble of
# posterior medians for all angular parameters (excluding the sampler's
# internal sigma_m scaling column).
extract_post_params_ang_star <- function(dep_type, dep_level, gauge, data_num) {
  params <- qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/%s/%s_%s_%s.qs",
                          dep_type, gauge, dep_level, data_num))$samples |>
    as_tibble()

  params |>
    select(-sigma_m) |>  # drop sampler bookkeeping column
    apply(MARGIN = 2, FUN = median) |>
    t() |>
    as_tibble() |>
    mutate(dataset = data_num)
}

# create_tib_post_params_ang_star --------------------------------------------
# Loops over all 200 datasets, stacks the one-row summaries, and saves.
create_tib_post_params_ang_star <- function(dep_type, dep_level, gauge) {
  tib_post_params <- sapply(
    1:200,
    function(x) extract_post_params_ang_star(dep_type, dep_level, gauge, x),
    simplify = FALSE
  ) |>
    bind_rows()

  filepath <- sprintf("fits_and_weights/post_params_joint/%s_%s_%s_ang_star.qs",
                      dep_type, dep_level, gauge)
  qsave(tib_post_params, filepath)
  print("Posterior medians of star-shaped angular parameters have been successfully saved")
}

create_tib_post_params_ang_star(dep_type = dep_type, dep_level = dep_level, gauge = gauge)
