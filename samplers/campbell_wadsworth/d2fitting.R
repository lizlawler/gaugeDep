# =============================================================================
# Fit the Campbell-Wadsworth (2023) piecewise-linear gauge model to the
# simulation-study datasets and compute exceedance-probability predictions for
# three prediction boxes. The radial model (fit.pwlin.2d) and angular model
# (fit.pwlin.2d with fW.fit = TRUE) are fit separately; predictions use an
# adaptive-k importance-sampling scheme to reach each box far in the tail.
#
# Called by: shell_scripts/run_all_mcmc_AD_wc.sh (or similar)
# Inputs:    data/{dep_type}/{dep_level}_{i}.json
# Outputs:   samplers/campbell_wadsworth/mle_and_preds/{dep_type}/{dep_level}_{i}.qs
#            (each file holds the preds tibble, radial_mle, and angular_mle)
#
# Command-line args:
#   1. dep_type   -- dependence structure ("gauss", "logistic", "husler_reiss")
#   2. dep_level  -- dependence strength ("low", "mid", "high")
#   3. data_start -- first dataset index
#   4. data_end   -- last dataset index
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
dep_type   <- args[1]
dep_level  <- args[2]
data_start <- args[3]
data_end   <- args[4]

library(geometry)
library(geometricMVE)
library(evd)
library(mvtnorm)
library(tidyr)
library(dplyr)
library(PWLExtremes)
library(qs)
source("samplers/campbell_wadsworth/mod_sim2d.R")

# Reference angles (knots) for the piecewise-linear gauge.
par.locs <- seq(0, 1, length.out = 11)
gfun <- function(w, par, locs = par.locs) gfun.2d(x = w, par = par, ref.angles = locs)

# Choose the importance-sampling extrapolation factor k for a box. The boxes
# sit far in the tail, so sampling at k = 1 rarely reaches them; scaling the
# threshold by k shifts the sampler's mass toward the box. Returns the smallest
# k (floored at 1) that covers a fine pseudo-grid tiling the box.
find_k <- function(box, ctau_val, radial_mle) {
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE        ~ c(2, 4)
  )

  pseudo_pred <- expand_grid(x1_pseudo = seq(dim1[1], dim1[2], length.out = 15),
                             x2_pseudo = seq(dim2[1], dim2[2], length.out = 15)) |>
    mutate(r_pseudo  = x1_pseudo + x2_pseudo,
           w1_pseudo = x1_pseudo / r_pseudo,
           w2_pseudo = x2_pseudo / r_pseudo)

  gw_pseudo <- sapply(pseudo_pred$w1_pseudo, function(x) gfun(x, radial_mle))
  poss_k    <- pseudo_pred$r_pseudo * gw_pseudo / ctau_val
  return(max(round(min(poss_k), 1) - 0.1, 1))
}

# Estimate P(X in box) by importance sampling from the fitted joint model.
# k > 1 returns points and importance weights; k = 1 returns points only.
# The estimate rescales by the mean weight and the exceedance rate to convert
# the conditional (above-threshold) probability to an unconditional one.
make_preds <- function(n = 50000, k = 1, box, radial_mle, angular_mle, par.locs, r, w, rexc) {
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE        ~ c(2, 4)
  )

  if (k > 1) {
    sim_list <- sim.2d.joint.mod(nsim = n, k.vals = k,
                                 gfun = gfun, par = radial_mle, fW.par = angular_mle,
                                 par.locs = par.locs, r = r, w = w)[[1]]
    xstar    <- sim_list$xstar
    r_over_k <- mean(sim_list$iw)
  } else {
    xstar <- sim.2d.joint.mod(nsim = n, k.vals = 1,
                              gfun = gfun, par = radial_mle, fW.par = angular_mle,
                              par.locs = par.locs, r = r, w = w)[[1]]
    r_over_k <- 1
  }

  prob.est <- mean(xstar[, 1] > dim1[1] & xstar[, 1] < dim1[2] &
                   xstar[, 2] > dim2[1] & xstar[, 2] < dim2[2]) *
    r_over_k * length(rexc) / length(r)
  return(prob.est)
}

for (i in data_start:data_end) {
  data <- RcppSimdJson::fload(sprintf("data/%s/%s_%s.json", dep_type, dep_level, i))

  r <- data$R
  w <- data$W
  x <- cbind(r * w, r * (1 - w))

  # Estimate the 95% radial threshold r0(w), then keep only the exceedances.
  qr     <- radial.quants.L1.KDE.2d(r, w, tau = 0.95, bww = 0.05, bwr = 0.05)
  r0w    <- qr$r0w
  wpts   <- qr$wpts
  excind <- r > r0w
  rexc   <- r[excind]
  wexc   <- w[excind]
  r0w    <- r0w[excind]

  # Fit the radial gauge (bounded) and the angular density separately.
  model.fit.R.bounded <- fit.pwlin.2d(r = rexc, r0w = r0w, w = wexc, locs = par.locs,
                                      pen.const = 1, method = "BFGS", bound.fit = TRUE)
  model.fit.W         <- fit.pwlin.2d(r = rexc, r0w = r0w, w = wexc, locs = par.locs,
                                      pen.const = NULL, fW.fit = TRUE, method = "BFGS")

  # Tail scaling constant from the fitted gauge, then predict each box.
  gw_fitted   <- sapply(w, function(x) gfun(x, model.fit.R.bounded$mle))
  ctau_fitted <- quantile(gw_fitted * r, 0.95)
  boxes       <- c("b1", "b2", "b3")
  preds_by_box <- sapply(boxes, function(x) {
    k    <- find_k(box = x, ctau_val = ctau_fitted, radial_mle = model.fit.R.bounded$mle)
    pred <- make_preds(n = 50000, k = k, box = x,
                       radial_mle = model.fit.R.bounded$mle, angular_mle = model.fit.W$fW.mle,
                       par.locs = par.locs, r = r, w = w, rexc = rexc)
    return(tibble(k = k, box = x, pred = pred, dataset = i))
  }, simplify = FALSE) |> bind_rows()

  results <- list(preds = preds_by_box,
                  radial_mle = model.fit.R.bounded$mle,
                  angular_mle = model.fit.W$fW.mle)
  qsave(x = results, file = sprintf("samplers/campbell_wadsworth/mle_and_preds/%s/%s_%s.qs",
                                    dep_type, dep_level, i))
  print(paste0("Successfully saved Campbell-Wadsworth fit for dataset number: ", i))
}
