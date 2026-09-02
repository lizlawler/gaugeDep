# =============================================================================
# Fit the star-shaped angular density on real fire-weather data at a single
# station. Runs three chains in parallel via furrr so convergence can be
# assessed with the posterior package.
#
# Called by: shell_scripts/local_machine/real_data_all_mcmc.sh (or similar)
# Inputs:    data/raw/{data_type}_expo.qs
# Outputs:   samplers/rcpp/ang_star_mcmc_fits/real_data/{data_type}_{gauge}.qs
#
# Command-line args:
#   1. gauge -- gauge function ("gauss", "logistic", "inv_log", etc.)
#
# Note: data_type is hardcoded to "friendmtn" below; change it to fit the
#       other station.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
gauge <- args[1]

library(furrr)
library(dplyr)
library(tidyr)
library(posterior)
library(qs)

plan(multisession, workers = 3)

data_type <- "friendmtn"  # station identifier; update for other stations

# Run a single angular MCMC chain. Wrapped for furrr; dependencies are
# re-loaded inside since each worker is a fresh R session.
run_chain <- function(chain_id, gauge, data_type) {
  library(gaugeDependence)
  library(qs)

  if (gauge != "dirichlet") {
    starting_vals <- runif(1)
  } else {
    starting_vals <- c(abs(rt(1, 4)) * 4, abs(rt(1, 4)) * 2)
  }

  fire_data_list <- qread(sprintf("data/raw/%s_expo.qs", data_type))
  w <- fire_data_list$W

  # More iterations than the simulation study (50k vs 15k) for the single
  # real dataset.
  results <- angular_mcmc(
    angles         = w,
    dim            = 2,
    starting_theta = starting_vals,
    gauge_type     = gauge,
    n_updates      = 50000,
    update_freq    = 250,
    adapt_cov      = TRUE
  )

  samples <- as.data.frame(results$samples) |>
    mutate(.chain = chain_id, .iteration = row_number(), .draw = NA_integer_)

  list(samples = samples, acceptance_rate = results$acc_prob)
}

all_chains <- future_map(
  1:3,
  ~ run_chain(.x, gauge = gauge, data_type = data_type),
  .options = furrr_options(seed = TRUE)
)

# Stack chains, assign a global draw index, and coerce to a draws_df.
combined_draws <- bind_rows(lapply(all_chains, function(x) x$samples)) |>
  group_by(.chain) |>
  mutate(.draw = row_number()) |>
  ungroup() |>
  as_draws_df()

print(sprintf("Now saving object to disk..."))
qsave(
  list(draws        = combined_draws,
       accept_rates = lapply(all_chains, function(x) x$acceptance_rate)),
  sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/%s_%s.qs", data_type, gauge)
)
print(sprintf("Successfully saved star-shaped MCMC fit on %s for %s gauge", data_type, gauge))
