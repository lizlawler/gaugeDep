library(tidyverse)
source("gauge_functions_wrt_w.R")

radial_cens_loglik <- function(radius, pars, threshold, w) {
  alpha <- pars[1]
  dep <- pars[2]
  beta <- gauss_gauge(w, dep)
  loglik <- 0.0
  for(i in seq_along(radius)) {
    if(radius[i] < threshold[i]) {
      loglik <- loglik + pgamma(threshold[i], shape = alpha, rate = beta[i], log.p = TRUE)
    } else {
      loglik <- loglik + dgamma(radius[i], shape = alpha, rate = beta[i], log = TRUE)
    }
  }
  return(loglik)
}

prior_lpdf <- function(pars, hyper_params) {
  alpha_prior <- dgamma(pars[1], shape = hyper_params[1], rate = hyper_params[2], log = TRUE)
  dep_prior <- dunif(pars[2], 0, 1, log = TRUE)
  return(alpha_prior + dep_prior)
}

adaptive.metr <- function(z, starting.theta,
                          likelihood.fn,
                          prior.fn,
                          hyper.params,
                          n.updates,
                          prop.Sigma = NULL,
                          adapt.cov = FALSE,
                          return.prop.Sigma.trace = FALSE,
                          r.opt = .234,
                          c.0 = 10,
                          c.1 = .8,
                          K = 10, ...) {
  
  call <- match.call()
  
  eps <- .001
  k <- 3  # the iteration offset
  
  
  p <- length(starting.theta)
  
  # If the supplied proposal covariance matrix is either not given or invalid,
  # just use the identity.
  if ((is.null(prop.Sigma)) |
      (length(prop.Sigma) != p^2) |
      (class(try(chol(prop.Sigma), silent=TRUE)) == "try-error")) {
    prop.Sigma <- diag(p)
    cat("Invalid or missing proposal covariance matrix.  Using identity.\n")
  }
  
  # Initialize sigma.m to the rule of thumb
  sigma.m <- (2.4/p)^2
  r.hat <- 0
  
  # Initialize prop.C
  prop.C <- chol(prop.Sigma)
  
  # Set up and initialize trace objects
  trace <- matrix(0, n.updates, p)
  sigma.m.trace <- rep(0, n.updates)
  r.trace <- rep(0, n.updates)
  jump.trace <- rep(0, n.updates)
  
  trace[1, ] <- starting.theta
  sigma.m.trace[1] <- sigma.m
  
  if (return.prop.Sigma.trace) {
    prop.Sigma.trace <- array(0, dim=c(n.updates, p, p))
    prop.Sigma.trace[1, , ] <- prop.Sigma
  }
  
  # Initialize Metropolis
  theta <- starting.theta
  likelihood <- likelihood.fn(z, theta, ...)
  prior <- prior.fn(theta, hyper.params)
  
  
  #########################################################
  # Begin main loop
  for (i in 2:n.updates) {
    theta.star <- theta + sigma.m * drop(rnorm(p) %*% prop.C)
    prior.star <- prior.fn(theta.star, hyper.params)
    if (prior.star != -Inf) {
      likelihood.star <- likelihood.fn(z, theta.star, ...)
      
      metr.ratio <- exp(prior.star + likelihood.star -
                          prior - likelihood)
      if (metr.ratio > runif(1)) {
        theta <- theta.star
        prior <- prior.star
        likelihood <- likelihood.star
        jump.trace[i] <- 1
      }
    }
    
    ########################################################
    # Adapt via my method                                  #
    if ((i %% K) == 0) {
      gamma1 <- c.0 / ((i/K) + k)^(c.1)
      gamma2 <- 1 / ((i/K) + k)^(c.1)
      
      r.hat <- mean(jump.trace[(i - K + 1) : i])
      
      sigma.m <- exp(log(sigma.m) +
                       gamma1*(r.hat - r.opt))
      
      if (adapt.cov) {
        prop.Sigma <- prop.Sigma +
          gamma2*(cov(trace[((i - K + 1) : i), ]) - prop.Sigma)
        
        while(is(try(chol(prop.Sigma), silent=TRUE), "try-error")) {
          prop.Sigma <- prop.Sigma + eps*diag(p)
          cat("Oops. Proposal covariance matrix is now:\n")
          print(prop.Sigma)
        }
        prop.C <- chol(prop.Sigma)
      }
    }
    
    #                                                      #
    ########################################################
    
    
    
    # Update the trace objects
    trace[i, ] <- theta
    sigma.m.trace[i] <- sigma.m
    r.trace[i] <- r.hat
    if (return.prop.Sigma.trace) {
      prop.Sigma.trace[i, , ] <- prop.Sigma
    }
    
    # Echo every 100 iterations
    if ((i %% 100) == 0)
      cat("Finished", i, "out of", n.updates, "iterations.\n")
  }
  # End main loop
  #########################################################
  
  # Collect trace objects to return
  res <- list(call=call,
              trace=trace,
              sigma.m.trace=sigma.m.trace,
              r.trace=r.trace,
              acc.prob=mean(jump.trace))
  if (return.prop.Sigma.trace) {
    res <- c(res, prop.Sigma.trace = list(prop.Sigma.trace))
  }
  
  return(res)
  
}

read_and_run <- function(dep_type, dep_level, data_num) {
  data <- RcppSimdJson::fload(paste0("data/", dep_type, "/", dep_level, "_", data_num, ".json"))
  R <- data$R
  r0w <- data$r0_w
  W <- data$W
}

samples <- adaptive.metr(z = R, starting.theta = c(rgamma(1, 4, 2), runif(1)), 
                         likelihood.fn = radial_cens_loglik, 
                         prior.fn = prior_lpdf, hyper.params = c(4,2), n.updates = 5000, 
                         threshold = r0w, w = W)

theta_samples <- samples$trace
sigma_trace <- samples$sigma.m.trace
accept_trace <- samples$r.trace
plot(accept_trace, type = "l")
plot(sigma_trace, type = "l")
plot(theta_samples[1000:5000,2], type = "l")

