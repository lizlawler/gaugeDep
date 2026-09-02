# =============================================================================
# R implementations of the extended GPD (EGPD) family (G1-G4) and the
# underlying GPD helper functions (pdf, CDF, CCDF, inverse CDF). Used for
# marginal threshold modelling of the fire weather data (ERC, FWI) via Stan.
# Sourced by: posterior_marg_transform.R, pred_task_real_data.R,
#             mix_dens_transform_fwi.R
# =============================================================================

## GPD building blocks --------------------------------------------------------
## Standard generalised Pareto pdf/cdf/ccdf/quantile. All are written on the log
## scale (log1p / log-sum-exp style) for numerical stability in the tail, with a
## `log` flag controlling whether the log or natural scale is returned.
gpd_pdf <- function(x, sigma, xi, log = FALSE) {
  lpdf <- -log(sigma) - (1/xi + 1) * log1p(xi * x / sigma)
  if(log) {
    return(lpdf)
  } else {
    return(exp(lpdf)) 
  }
}

gpd_cdf <- function(x, sigma, xi, log = FALSE) {
  lcdf <- log1p(-exp(-1/xi * log1p(xi * x / sigma))) 
  if(log) {
    return(lcdf)
  } else {
    return(exp(lcdf))
  }
}

gpd_ccdf <- function(x, sigma, xi, log = FALSE) {
  lccdf <- -1/xi * log1p(xi * x/sigma)
  if(log) {
    return(lccdf)
  } else {
    return(exp(lccdf))
  }
}

gpd_icdf <- function(q, sigma, xi) {
  return(sigma / xi * ((1 - q)^(-xi) - 1))
}

## EGPD families -------------------------------------------------------------
## Extended GPD families G1-G4 (Naveau et al. 2016). Each applies a different
## carrier transform to the GPD cdf so the model has GPD tails but flexible
## behaviour near the lower end, avoiding a hard threshold choice.

# G1: power-law carrier, single shape kappa.
g1_pdf <- function(x, sigma, xi, kappa, log = FALSE) {
  lpdf <- log(kappa) + 
    (kappa - 1) * gpd_cdf(x, sigma, xi, TRUE) + 
    gpd_pdf(x, sigma, xi, TRUE)
  if(log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

g1_cdf <- function(x, sigma, xi, kappa, log = FALSE) {
  lcdf <- kappa * gpd_cdf(x, sigma, xi, TRUE)
  if(log) {
    return(lcdf)
  } else {
    return(exp(lcdf))
  }
}

g1_icdf <- function(q, sigma, xi, kappa) {
  icdf <- gpd_icdf(q ^ (1/kappa), sigma, xi)
  return(icdf)
}

# G2: two-component mixture of G1 carriers with weight p.
g2_pdf <- function(x, sigma, xi, p, kappa1, kappa2, log = FALSE) {
  lpdf <- gpd_pdf(x, sigma, xi, TRUE) + 
    (kappa1 - kappa2) * gpd_cdf(x, sigma, xi, TRUE) +
    log(
      exp(log(kappa1) + log(p)) +
        exp(log(kappa2) + log1p(-p))
    )
  if(log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

g2_cdf <- function(x, sigma, xi, p, kappa1, kappa2, log = FALSE) {
  lcdf <- log(
    exp(log(p) + g1_cdf(x, sigma, xi, kappa1, TRUE)) +
      exp(log1p(-p) + g1_cdf(x, sigma, xi, kappa2, TRUE))
    )
  if(log) {
    return(lcdf)
  } else {
    return(exp(lcdf))
  }
} 

# G3: Beta carrier with parameter delta.
g3_pdf <- function(x, sigma, xi, delta, log = FALSE) {
  inner_x <- exp(delta * gpd_ccdf(x, sigma, xi, TRUE))
  lpdf <- dbeta(inner_x, shape1 = 1/delta, shape2 = 2, log = TRUE) + 
    log(delta) + 
    (delta - 1) * gpd_ccdf(x, sigma, xi, TRUE) +
    gpd_pdf(x, sigma, xi, TRUE)
  if(log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

g3_cdf <- function(x, sigma, xi, delta, log = FALSE) {
  inner_x <- exp(delta * gpd_ccdf(x, sigma, xi, TRUE))
  lcdf <- pbeta(inner_x, shape1 = 1/delta, shape2 = 2, lower.tail = FALSE, log.p = TRUE)
  if(log) {
    return(lcdf)
  } else {
    return(exp(lcdf))
  }
}

# G4: G3 carrier raised to a further power kappa/2.
g4_pdf <- function(x, sigma, xi, delta, kappa, log = FALSE) {
  lpdf <- log(kappa) - log(2) + (kappa / 2 - 1) * g3_cdf(x, sigma, xi, delta, TRUE) + g3_pdf(x, sigma, xi, delta, TRUE)
  if(log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

g4_cdf <- function(x, sigma, xi, delta, kappa, log = FALSE) {
  lcdf <- kappa / 2 * g3_cdf(x, sigma, xi, delta, TRUE)
  if(log) {
    return(lcdf)
  } else {
    return(exp(lcdf))
  }
}
