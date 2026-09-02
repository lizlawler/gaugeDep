# =============================================================================
# Produce MCMC traceplots for the Stan marginal-transform fits, for visual
# convergence assessment. Reads the CmdStan CSV output for a given fit and
# writes per-parameter traceplots via MCMCvis::MCMCtrace().
#
# Inputs:  samplers/stan/marg_transform/csv_fits/  (CmdStan CSV output)
# Outputs: traceplot PDF(s)
#
# Note: set the `pattern` below to match the fit you want to inspect (the
#       example targets the "thomescreek_g2" marginal fit).
# =============================================================================

library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)

# Locate the CmdStan CSV files for the fit of interest and load them as an
# mcmc.list for MCMCvis.
files <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                    pattern = "thomescreek_g2", full.names = TRUE)

fit <- as_mcmc.list(as_cmdstan_fit(files))

# One traceplot per parameter, opened in the PDF viewer.
MCMCtrace(fit,
          ind      = TRUE,
          open_pdf = TRUE)
