# =============================================================================
# Computes BMA model weights (stacking, pseudo-BMA, pseudo-BMA+) for the joint
# model (radial + angular components combined) across all 200 simulation study
# datasets. The joint log-likelihood for each dataset is assembled from the
# pre-computed per-component pointwise log-likelihood matrices, then passed to
# loo::stacking_weights() and loo::pseudobma_weights().
#
# Called by: shell_scripts/run_joint_model_wts.sh
# Inputs:    samplers/rcpp/radial_mcmc_fits/{dep_type}/pw_loglik/...qs
#            samplers/rcpp/ang_star_mcmc_fits/{dep_type}/pw_loglik/...qs  (star)
#            samplers/nimble/ang_mix_mcmc_fits/{dep_type}/pw_loglik/...qs (mix)
#            data/{dep_type}/{dep_level}_{i}.json                         (for idx)
# Outputs:   fits_and_weights/wts_joint_model/{dep_type}_{likelihood}_{angle_dens}_{dep_level}.qs
#            samplers/joint_loglik/{dep_type}/{angle_dens}/...qs           (cached)
#
# Command-line args:
#   1. dep_type   -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level  -- dependence strength ("low", "mid", "high")
#   3. likelihood -- likelihood type ("trunc" or "cens")
#   4. angle_dens -- angular density model ("star" or "mix")
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type   <- args[1]
dep_level  <- args[2]
likelihood <- args[3]
angle_dens <- args[4]

library(qs)
library(posterior)
library(loo)
library(dplyr)
library(tidyr)

options(mc.cores = parallel::detectCores())

# create_joint_loglik --------------------------------------------------------
# Assembles the joint (radial + angular) pointwise log-likelihood matrix for a
# single dataset by summing the component matrices. Caches the result to disk
# to avoid redundant computation across weight-extraction runs.
create_joint_loglik <- function(dep_type, dep_level, gauge, data_num,
                                likelihood = "cens", angle_dens = "star") {
  loglik_file <- sprintf("samplers/joint_loglik/%s/%s/%s_%s_%s_%s.qs",
                         dep_type, angle_dens, gauge, likelihood, dep_level, data_num)

  # Return cached file if it exists
  if (file.exists(loglik_file)) return(qread(loglik_file))

  temp_radial <- qread(sprintf(
    "samplers/rcpp/radial_mcmc_fits/%s/pw_loglik/%s_%s_%s_%s.qs",
    dep_type, gauge, likelihood, dep_level, data_num
  ))

  temp_angular <- if (angle_dens == "star") {
    qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/%s/pw_loglik/%s_%s_%s.qs",
                  dep_type, gauge, dep_level, data_num))
  } else {
    qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/pw_loglik/%s_%s.qs",
                  dep_type, dep_level, data_num))
  }

  # For truncated likelihood, restrict angular loglik to exceedance indices only
  if (likelihood == "trunc") {
    idx <- RcppSimdJson::fload(
      sprintf("data/%s/%s_%s.json", dep_type, dep_level, data_num)
    )$idx
    temp_angular <- temp_angular[, idx]
  }

  temp_joint <- temp_radial + temp_angular
  qsave(temp_joint, loglik_file)
  temp_joint
}

# extract_lpd_pt -------------------------------------------------------------
# Runs LOO-CV on the joint loglik matrix for one dataset and gauge, returning
# the vector of observation-level expected log predictive densities (ELPD).
extract_lpd_pt <- function(dep_type, dep_level, gauge, data_num, likelihood, angle_dens) {
  temp     <- create_joint_loglik(dep_type, dep_level, gauge, data_num, likelihood, angle_dens)
  loo_temp <- loo(temp)
  loo_temp$pointwise[, "elpd_loo"]
}

# create_lpd_list ------------------------------------------------------------
# Returns a named list of ELPD vectors (one per gauge) for a single dataset,
# ready to pass to stacking_weights() / pseudobma_weights().
create_lpd_list <- function(dep_type, dep_level, data_num, likelihood, angle_dens) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  setNames(
    sapply(gauge_library,
           function(x) extract_lpd_pt(dep_type, dep_level, x, data_num, likelihood, angle_dens)),
    gauge_library
  )
}

# model_weights --------------------------------------------------------------
# Computes all three BMA weight variants for a single dataset.
model_weights <- function(dep_type, dep_level, data_num, likelihood, angle_dens) {
  temp <- create_lpd_list(dep_type, dep_level, data_num, likelihood, angle_dens)
  print(paste0("Weights extracted for dataset: ", data_num))
  list(
    stacking        = stacking_weights(temp),
    pseudobma_boot  = pseudobma_weights(temp),
    pseudobma_noboot = pseudobma_weights(temp, BB = FALSE)
  )
}

# make_wts_df ----------------------------------------------------------------
# Iterates over all 200 datasets, collects BMA weights, reshapes into a tidy
# tibble, and saves. Short-circuits if the output file already exists.
make_wts_df <- function(dep_type, dep_level, likelihood, angle_dens) {
  wts_file <- sprintf("fits_and_weights/wts_joint_model/%s_%s_%s_%s.qs",
                      dep_type, likelihood, angle_dens, dep_level)

  if (file.exists(wts_file)) return(qread(wts_file))

  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")

  wts <- lapply(1:200, function(x) {
    model_weights(dep_type, dep_level, x, likelihood, angle_dens)
  }) |>
    bind_rows() |>
    mutate(method   = rep(gauge_library, 200),
           stacking = as.numeric(stacking),
           pseudobma_boot   = as.numeric(pseudobma_boot),
           pseudobma_noboot = as.numeric(pseudobma_noboot),
           dataset  = rep(1:200, each = 6))

  qsave(wts, wts_file)
  print(sprintf("Model weights for %s, %s, %s, %s have been created and saved to disk",
                dep_type, dep_level, likelihood, angle_dens))
}

make_wts_df(dep_type = dep_type, dep_level = dep_level,
            likelihood = likelihood, angle_dens = angle_dens)
