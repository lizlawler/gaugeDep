#########################
##### GPD functions #####
#########################
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
  licdf <- log(sigma) - log(xi) + log(exp(-xi * log1p(-q)) - 1)
  return(exp(licdf))
}

##########################
##### EGPD functions #####
##########################

## G1
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

## G2
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

## G3
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

## G4
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
