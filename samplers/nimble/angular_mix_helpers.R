# =============================================================================
# Shared helper functions for computing pointwise log-likelihoods under the
# angular mixture (stick-breaking / Dirichlet process mixture of Betas) model
# fit via NIMBLE. Sourced by:
#   - ang_mix_loglik_calc.R          (simulation study)
#   - ang_mix_loglik_calc_real_data.R (real data analysis)
# =============================================================================

# mix_lpdf -------------------------------------------------------------------
# Evaluates the mixture-of-Betas density at a vector of angles using a single
# set of mixture weights and Beta component parameters.
#
# Args:
#   w    - numeric vector of observed angles in (0, 1)
#   wts  - numeric vector of mixture weights (length L, summing to 1)
#   alpha - numeric vector of Beta shape1 parameters (length L)
#   beta  - numeric vector of Beta shape2 parameters (length L)
#
# Returns: numeric vector of mixture density values, same length as w
mix_lpdf <- function(w, wts, alpha, beta) {
  n    <- length(wts)
  dens <- 0.0
  for (i in 1:n) {
    dens <- dens + wts[i] * dbeta(w, alpha[i], beta[i])
  }
  return(dens)
}

# angular_loglik -------------------------------------------------------------
# Computes a posterior pointwise log-likelihood matrix for the angular mixture
# model. For each MCMC iteration, evaluates the mixture density at all observed
# angles. The resulting (n_iter x n_obs) matrix is passed to loo::loo() or
# loo::stacking_weights() for model comparison and BMA weight computation.
#
# Args:
#   angles           - numeric vector of observed angles in (0, 1)
#   posterior_params - matrix of posterior samples (rows = iterations),
#                      with named columns from NIMBLE output including
#                      "probs[k]", "alphastar[k]", and "betastar[k]"
#
# Returns: numeric matrix of dimension (n_iter x n_obs)
angular_loglik <- function(angles, posterior_params) {
  probs <- posterior_params[, grepl(pattern = "probs\\[",     colnames(posterior_params))]
  alpha <- posterior_params[, grepl(pattern = "alphastar\\[", colnames(posterior_params))]
  beta  <- posterior_params[, grepl(pattern = "betastar\\[",  colnames(posterior_params))]

  n_iter <- nrow(probs)
  n_obs  <- length(angles)

  pw_loglik <- matrix(NA, nrow = n_iter, ncol = n_obs)
  for (i in 1:n_iter) {
    pw_loglik[i, ] <- mix_lpdf(angles, probs[i, ], alpha[i, ], beta[i, ])
  }
  return(pw_loglik)
}
