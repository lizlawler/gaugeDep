# =============================================================================
# Fit the radial component on real fire-weather data (ERC/FWI at a single
# station). Runs three chains in parallel via furrr so that convergence can be
# assessed with the posterior package.
#
# Called by: shell_scripts/local_machine/real_data_all_mcmc.sh (or similar)
# Inputs:    data/raw/{data_type}_expo.qs
# Outputs:   samplers/rcpp/radial_mcmc_fits/real_data/
#              {data_type}_{gauge}_{likelihood}.qs
#
# Command-line args:
#   1. gauge      -- gauge function ("gauss", "logistic", "inv_log", etc.)
#   2. likelihood -- likelihood type ("trunc" or "cens")
#
# Note: data_type is hardcoded to "friendmtn" below; change it (or promote it
#       to a command-line argument) to fit the other station.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
gauge      <- args[1]
likelihood <- args[2]

library(furrr)
library(dplyr)
library(tidyr)
library(posterior)
library(qs)

plan(multisession, workers = 3)

data_type <- "friendmtn"  # station identifier; update for other stations

# Run a single MCMC chain for the radial model. Wrapped in a function so furrr
# can dispatch the three chains to separate workers; because those workers are
# fresh R sessions, the package dependencies are re-loaded inside the function.
run_chain <- function(chain_id, gauge, likelihood, data_type) {
  library(gaugeDependence)
  library(qs)

  if (gauge != "dirichlet") {
    starting_vals <- c(rgamma(1, 4, 2), runif(1))
  } else {
    starting_vals <- c(rgamma(1, 4, 2),
                       abs(rt(1, 4)) * 4,
                       abs(rt(1, 4)) * 2)
  }

  fire_data_list <- qread(sprintf("data/raw/%s_expo.qs", data_type))
  w   <- fire_data_list$W
  r   <- fire_data_list$R
  r0w <- fire_data_list$r0_w
  idx <- fire_data_list$idx

  if (likelihood == "trunc") {
    w   <- w[idx]
    r   <- r[idx]
    r0w <- r0w[idx]
  }

  # More iterations than the simulation study (50k vs 15k) since the single
  # real dataset is cheap to fit and benefits from the longer run.
  results <- radial_adaptive_mh(
    radii           = r,
    r0w             = r0w,
    angles          = w,
    starting_theta  = starting_vals,
    likelihood_type = likelihood,
    gauge_type      = gauge,
    n_updates       = 50000,
    update_freq     = 250,
    adapt_cov       = TRUE
  )

  # Tag each draw with chain / iteration metadata for the posterior package.
  samples <- as.data.frame(results$samples) |>
    mutate(.chain = chain_id, .iteration = row_number(), .draw = NA_integer_)

  list(samples = samples, acceptance_rate = results$acc_prob)
}

all_chains <- future_map(
  1:3,
  ~ run_chain(.x, gauge = gauge, likelihood = likelihood, data_type = data_type),
  .options = furrr_options(seed = TRUE)
)

# Stack the three chains and assign a global draw index, then coerce to a
# draws_df. (as_draws_df() acts on the piped object and takes no argument.)
combined_draws <- bind_rows(lapply(all_chains, function(x) x$samples)) |>
  group_by(.chain) |>
  mutate(.draw = row_number()) |>
  ungroup() |>
  as_draws_df()

print(sprintf("Now saving object to disk..."))
qsave(
  list(draws        = combined_draws,
       accept_rates = lapply(all_chains, function(x) x$acceptance_rate)),
  sprintf("samplers/rcpp/radial_mcmc_fits/real_data/%s_%s_%s.qs",
          data_type, gauge, likelihood)
)
print(sprintf(
  "Successfully saved radial MCMC fit on fire data for %s, %s gauge, %s likelihood",
  data_type, gauge, likelihood
))
