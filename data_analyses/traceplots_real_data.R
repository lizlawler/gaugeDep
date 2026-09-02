# =============================================================================
# Interactive convergence assessment for the real-data MCMC fits. Traceplots
# and posterior densities are inspected on screen for each gauge/likelihood to
# decide which chain to discard; the retained two chains are then written out
# as the *_2chains.qs objects that the log-likelihood and prediction scripts
# consume. Run this script interactively, one block at a time.
#
# Inputs:    samplers/rcpp/radial_mcmc_fits/real_data/{data_type}_*.qs
#            samplers/rcpp/ang_star_mcmc_fits/real_data/{data_type}_*.qs
#            samplers/nimble/ang_mix_mcmc_fits/real_data/{data_type}.qs
# Outputs:   samplers/rcpp/radial_mcmc_fits/real_data/{data_type}_{gauge}_{lhood}_2chains.qs
#            samplers/rcpp/ang_star_mcmc_fits/real_data/{data_type}_{gauge}_2chains.qs
#            samplers/nimble/ang_mix_mcmc_fits/real_data/{data_type}_2chains.qs
# =============================================================================

library(qs)
library(tidyverse)
library(posterior)
library(grafify)
library(patchwork)
library(bayesplot)
library(coda)

data_type <- "friendmtn"

qs_files <- list.files(path = "samplers/rcpp/radial_mcmc_fits/real_data/", pattern = data_type, full.names = TRUE)
qs_files <- qs_files[!grepl("2chains", qs_files)]
object_name <- str_remove(str_remove(basename(qs_files), sprintf("%s_", data_type)), ".qs")

# Load each fit's draws into the global environment under a name derived from
# its filename (e.g. "gauss_trunc"), so the plotting helpers below can look
# them up with get().
read_files <- function(filepath, object) {
  temp <- qread(filepath)$draws
  assign(object, temp, envir = .GlobalEnv)
  rm(temp)
  gc()
}

for (i in seq_along(qs_files)) {
  read_files(qs_files[i], object_name[i])
}

gg_trace_plot <- function(fit_mcmc, gauge) {
  if (gauge == "dirichlet") {
    p1 <- fit_mcmc |> ggplot(aes(x = .iteration, y = theta1, color = as.factor(.chain))) +
      geom_line(alpha = 0.75) +
      scale_color_grafify(palette = "r4", ColSeq = FALSE) +
      theme_classic()
    p2 <- fit_mcmc |> ggplot(aes(x = .iteration, y = theta2, color = as.factor(.chain))) +
      geom_line(alpha = 0.75) +
      scale_color_grafify(palette = "r4", ColSeq = FALSE) +
      theme_classic()
    return(list(t1 = p1, t2 = p2))
  } else {
    fit_mcmc |> ggplot(aes(x = .iteration, y = dep, color = as.factor(.chain))) +
      geom_line(alpha = 0.75) +
      scale_color_grafify(palette = "r4", ColSeq = FALSE) +
      theme_classic()
  }
}

gg_dens_plot <- function(fit_mcmc, gauge) {
  if (gauge == "dirichlet") {
    p1 <- fit_mcmc |> ggplot(aes(x = theta1, color = as.factor(.chain))) +
      geom_density() +
      scale_color_grafify(palette = "r4", ColSeq = FALSE, guide = "none") +
      theme_classic()
    p2 <- fit_mcmc |> ggplot(aes(x = theta2, color = as.factor(.chain))) +
      geom_density() +
      scale_color_grafify(palette = "r4", ColSeq = FALSE, guide = "none") +
      theme_classic()
    return(list(d1 = p1, d2 = p2))
  } else {
    fit_mcmc |> ggplot(aes(x = dep, color = as.factor(.chain))) +
      geom_density() +
      scale_color_grafify(palette = "r4", ColSeq = FALSE, guide = "none") +
      theme_classic()
  }
}

# Show density + trace for all iterations (top row) against a thinned,
# post-warmup subset with optional chains dropped (bottom row), so the effect
# of removing a chain is visible before committing to it.
trace_and_dens_plots <- function(gauge, radial = FALSE, likelihood = NULL, idx_from = 30001, idx_to = 50000, thin = 10, chain_remove = NULL) {
  iter_idx <- seq.int(idx_from, idx_to, by = thin)

  if (radial) {
    fit_obj <- get(sprintf("%s_%s", gauge, likelihood))
  } else {
    fit_obj <- get(gauge)
  }

  all_iter_trace <- gg_trace_plot(fit_obj, gauge)
  all_iter_dens <- gg_dens_plot(fit_obj, gauge)
  if (!is.null(chain_remove)) {
    subset_trace <- gg_trace_plot((fit_obj |> filter(.iteration %in% iter_idx, !.chain %in% chain_remove)), gauge)
    subset_dens <- gg_dens_plot((fit_obj |> filter(.iteration %in% iter_idx, !.chain %in% chain_remove)), gauge)
  } else {
    subset_trace <- gg_trace_plot((fit_obj |> filter(.iteration %in% iter_idx)), gauge)
    subset_dens <- gg_dens_plot((fit_obj |> filter(.iteration %in% iter_idx)), gauge)
  }

  if (gauge == "dirichlet") {
    p1 <- all_iter_dens$d1 + all_iter_trace$t1
    p2 <- subset_dens$d1 + subset_trace$t1
    p3 <- all_iter_dens$d2 + all_iter_trace$t2
    p4 <- subset_dens$d2 + subset_trace$t2
    return(p1 / p2 / p3 / p4)
  } else {
    return((all_iter_dens + all_iter_trace) / (subset_dens + subset_trace))
  }
}

trace_and_dens_plots("gauss", radial = TRUE, likelihood = "trunc")

## Friend Mtn chain selection (removal)
## dirichlet: cens, 2; trunc, 2
## asym_log: cens, 1; trunc, 2
## logistic: cens, 3; trunc, 3
## rectangular: cens, 2; trunc, 3
## inv_log: cens, 2; trunc, 3
## gauss: cens, 3; trunc, 3

## Redstone chain selection (removal)
## asym_log, trunc and cens, chain 1; logistic, cens, chain 2
## since we only have two chains for the above, pick the best two chains from all others to simplify weighting schemas
## gauss cens and trunc: 3,
## inv_log cens:1, trunc:3
## rectangular trunc: 1, cens: 2
## dirichlet cens: 1, trunc: 1
## logistic trunc: 1

# Write out the retained chains as a plain parameter matrix for downstream use.
subset_iters_remove_chain_radial <- function(fit_obj, chain_remove) {
  sub_idx <- seq.int(30001, 50000, by = 10)
  temp <- fit_obj |> filter(.iteration %in% sub_idx, !.chain %in% chain_remove)
  if (grepl("dirichlet", deparse(substitute(fit_obj)))) {
    temp <- temp |>
      select(alpha, theta1, theta2) |>
      as.matrix()
  } else {
    temp <- temp |>
      select(alpha, dep) |>
      as.matrix()
  }
  qsave(temp, sprintf("samplers/rcpp/radial_mcmc_fits/real_data/%s_%s_2chains.qs", data_type, deparse(substitute(fit_obj))))
}

subset_iters_remove_chain_radial(gauss_trunc, 3)

## angular densities
qs_files <- list.files(path = "samplers/rcpp/ang_star_mcmc_fits/real_data/", pattern = data_type, full.names = TRUE)
qs_files <- qs_files[!grepl("2chains", qs_files)]
object_name <- str_remove(str_remove(basename(qs_files), sprintf("%s_", data_type)), ".qs")

for (i in seq_along(qs_files)) {
  read_files(qs_files[i], object_name[i])
}

subset_iters_remove_chain_angular <- function(fit_obj, chain_remove, idx_from = 30001, idx_to = 50000, thin = 10) {
  sub_idx <- seq.int(idx_from, idx_to, by = thin)
  temp <- fit_obj |> filter(.iteration %in% sub_idx, !.chain %in% chain_remove)
  if (grepl("dirichlet", deparse(substitute(fit_obj)))) {
    temp <- temp |>
      select(theta1, theta2) |>
      as.matrix()
  } else {
    temp <- temp |>
      select(dep) |>
      as.matrix()
  }
  qsave(temp, sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/%s_%s_2chains.qs", data_type, deparse(substitute(fit_obj))))
}

trace_and_dens_plots("dirichlet", chain_remove = 1)

subset_iters_remove_chain_angular(dirichlet, chain_remove = 1)

mix_mcmc <- qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/%s.qs", data_type))

# Convert nimble object to coda::mcmc.list if needed
mcmc_list <- as.mcmc.list(lapply(mix_mcmc, mcmc))

# Plot all alphastar[i] traces in one faceted plot
mcmc_trace(mcmc_list[c(1, 2)], regex_pars = "probs\\[")

mix_2chains <- cbind(rbind(mix_mcmc[[1]], mix_mcmc[[2]]), chain_id = rep(1:2, c(2000, 2000)))
qsave(mix_2chains, file = sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/%s_2chains.qs", data_type))
