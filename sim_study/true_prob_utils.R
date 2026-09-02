# =============================================================================
# Utility functions for computing true bivariate exceedance probabilities
# under the simulation study's data-generating models. Sourced by all
# sim_study/ scripts that need ground-truth probabilities for evaluation.
#
# Depends on: mvtnorm, evd
# =============================================================================

# true_gauss_prob ------------------------------------------------------------
# Computes P(X1 in dim1, X2 in dim2) where (X1, X2) follow a bivariate
# Gaussian copula with standard exponential margins, using the probability
# integral transform to Gaussian scale.
#
# Args:
#   dim1 - numeric vector of length 2: c(lower, upper) bound for X1
#   dim2 - numeric vector of length 2: c(lower, upper) bound for X2
#   dep  - scalar correlation parameter in (-1, 1)
# Returns: scalar probability
true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(mvtnorm::pmvnorm(lower = c(dim1_star[1], dim2_star[1]),
                          upper = c(dim1_star[2], dim2_star[2]),
                          corr  = corr_matrix)[1])
}

# true_bvevd_prob ------------------------------------------------------------
# Computes P(X1 in dim1, X2 in dim2) where (X1, X2) follow a bivariate
# extreme value distribution (logistic or Husler-Reiss) with standard
# exponential margins. Margins are transformed to standard Gumbel scale
# before applying the BVEVD CDF via inclusion-exclusion.
#
# Args:
#   dim1       - numeric vector of length 2: c(lower, upper) bound for X1
#   dim2       - numeric vector of length 2: c(lower, upper) bound for X2
#   dep        - scalar dependence parameter (see evd::pbvevd for interpretation)
#   model_type - character string for evd::pbvevd, e.g. "log" or "hr"
# Returns: scalar probability
true_bvevd_prob <- function(dim1, dim2, dep, model_type) {
  dim1_star <- evd::qgev(pexp(dim1), loc = 0, scale = 1, shape = 0)
  dim2_star <- evd::qgev(pexp(dim2), loc = 0, scale = 1, shape = 0)
  upper_right <- evd::pbvevd(q = c(dim1_star[2], dim2_star[2]), model = model_type, dep = dep)
  upper_left  <- evd::pbvevd(q = c(dim1_star[1], dim2_star[2]), model = model_type, dep = dep)
  lower_right <- evd::pbvevd(q = c(dim1_star[2], dim2_star[1]), model = model_type, dep = dep)
  lower_left  <- evd::pbvevd(q = c(dim1_star[1], dim2_star[1]), model = model_type, dep = dep)
  return(upper_right - upper_left - lower_right + lower_left)
}

# lookup_true_dep ------------------------------------------------------------
# Centralises the dependence parameter lookup that was previously copy-pasted
# across every sim_study script. Returns a list with elements `dep` (the
# scalar dependence parameter) and `model_type` (for use with true_bvevd_prob).
#
# Args:
#   dep_type  - character: "gauss", "logistic", or "husler_reiss"
#   dep_level - character: "low", "mid", "high", or _wc variants where defined
# Returns: named list with elements `dep` and `model_type`
lookup_true_dep <- function(dep_type, dep_level) {
  if (dep_type == "gauss") {
    levels_list <- list(low = 0.1, mid = 0.5, high = 0.9, high_wc = 0.8)
    return(list(dep = as.numeric(levels_list[[dep_level]]), model_type = "gauss"))
  } else if (dep_type == "logistic") {
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1, low_wc = 0.8, mid_wc = 0.4)
    return(list(dep = as.numeric(levels_list[[dep_level]]), model_type = "log"))
  } else {  # husler_reiss
    levels_list <- list(low = 0.1, mid = 1, high = 3)
    return(list(dep = as.numeric(levels_list[[dep_level]]), model_type = "hr"))
  }
}

# get_true_prob --------------------------------------------------------------
# Convenience wrapper combining lookup_true_dep with the appropriate
# probability function.
#
# Args:
#   dep_type  - character: "gauss", "logistic", or "husler_reiss"
#   dep_level - character: "low", "mid", or "high"
#   dim1      - numeric vector of length 2: c(lower, upper) bound for X1
#   dim2      - numeric vector of length 2: c(lower, upper) bound for X2
# Returns: scalar true probability
get_true_prob <- function(dep_type, dep_level, dim1, dim2) {
  params <- lookup_true_dep(dep_type, dep_level)
  if (dep_type == "gauss") {
    return(true_gauss_prob(dim1, dim2, params$dep))
  } else {
    return(true_bvevd_prob(dim1, dim2, params$dep, params$model_type))
  }
}
