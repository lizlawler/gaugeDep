# =============================================================================
# Extracts and summarises posterior parameters from the NIMBLE angular mixture
# (stick-breaking DP mixture of Betas) MCMC fits across all 200 simulation
# study datasets. Saves posterior means of mixture weights (probs) and Beta
# component parameters (alphastar, betastar) for downstream visualisation and
# prediction.
#
# Called by: shell_scripts/call_ang_mix_params.sh
# Inputs:    samplers/nimble/ang_mix_mcmc_fits/{dep_type}/{dep_level}_{i}.qs
# Outputs:   fits_and_weights/post_params_joint/{dep_type}_{dep_level}_ang_mix.qs
#
# Command-line args:
#   1. dep_type  -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level -- dependence strength ("low", "mid", "high")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type  <- args[1]
dep_level <- args[2]

library(qs)
library(dplyr)
library(tidyr)

# extract_post_params_ang_mix ------------------------------------------------
# Loads NIMBLE MCMC output for one dataset and returns a one-row tibble of
# posterior means for the mixture weights and Beta shape parameters, plus a
# dataset index column.
extract_post_params_ang_mix <- function(dep_type, dep_level, data_num) {
  params <- qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
                          dep_type, dep_level, data_num)) |>
    as_tibble() |>
    select(matches("probs|alphastar|betastar"))

  params |>
    colMeans() |>
    t() |>
    as_tibble() |>
    mutate(dataset = data_num)
}

# create_tib_post_params_ang_mix ---------------------------------------------
# Loops over all 200 datasets, stacks the one-row summaries into a single
# tibble, and saves to disk.
create_tib_post_params_ang_mix <- function(dep_type, dep_level) {
  tib_post_params <- sapply(
    1:200,
    function(x) extract_post_params_ang_mix(dep_type, dep_level, x),
    simplify = FALSE
  ) |>
    bind_rows()

  filepath <- sprintf("fits_and_weights/post_params_joint/%s_%s_ang_mix.qs",
                      dep_type, dep_level)
  qsave(tib_post_params, filepath)
  print("Posterior means of stick-breaking angular parameters have been successfully saved")
}

create_tib_post_params_ang_mix(dep_type = dep_type, dep_level = dep_level)
