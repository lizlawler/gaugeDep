sim.2d.joint.mod <- function(nsim, k.vals = 1, gfun, shape = 2, par, fW.par, par.locs, 
                             r, w, tau = 0.95, bww = 0.05, bwr = 0.05, marg = "pos") {
  if (marg == "Rd") {
    w.beta = (w + 2)/4
    beta.nll = function(pars) {
      -sum(dbeta(x = w.beta, shape1 = pars[1], shape2 = pars[2], 
                 log = T))
    }
  }
  else {
    beta.nll = function(pars) {
      -sum(dbeta(x = w, shape1 = pars[1], shape2 = pars[2], 
                 log = T))
    }
  }
  beta.mle = optim(par = rep(1, 2), fn = beta.nll)$par
  nthin = 2
  wstar = fW.mcmc.g.2d(niter = (nthin * nsim) + 1000, nburn = 1000, 
                       bpar1 = beta.mle[1], bpar2 = beta.mle[2], thin = nthin, 
                       g = gfun, par = fW.par, locs = par.locs)
  if (marg == "Rd") {
    wstar = wstar * 4 - 2
  }
  r.tau.wstar = sapply(wstar, function(wstar.i) {
    weightsw <- dnorm(w, mean = as.numeric(wstar.i), sd = bww)
    ccdf <- function(rc) {
      mean(weightsw * pnorm(rc, mean = r, sd = bwr))/mean(weightsw)
    }
    dummy <- function(rc) {
      ccdf(rc) - tau
    }
    ur <- uniroot(dummy, interval = c(0, 30))
    return(ur$root)
  })
  sims = lapply(k.vals, function(k) {
    if (k == 1) {
      rate <- sapply(wstar, gfun, par = par)
      rstar <- qgamma(1 - runif(nsim) * pgamma(k * r.tau.wstar, 
                                               shape = shape, rate = rate, lower.tail = F), 
                      shape = shape, rate = rate)
      if (marg == "pos") {
        xstar <- cbind(rstar * wstar, rstar * (1 - wstar))
      }
      else if (marg == "Rd") {
        xstar <- pol2cart.L1Rd(r = rstar, w = wstar)
      }
      return(xstar)
    }
    else if (k > 1) {
      iw <- iweights.2d.pwl(k = k, r0w = r.tau.wstar, w = wstar, 
                            gfun = gfun, shape = shape, par = par)
      star.ind <- sample(1:length(wstar), size = nsim, 
                         replace = T, prob = iw)
      wstar2 <- wstar[star.ind]
      r.tau.wstar2 <- c(k * r.tau.wstar)[star.ind]
      rate <- sapply(wstar2, gfun, par = par)
      rstar <- qgamma(1 - runif(nsim) * pgamma(r.tau.wstar2, 
                                               shape = shape, rate = rate, lower.tail = F), 
                      shape = shape, rate = rate)
      if (marg == "pos") {
        xstar <- cbind(rstar * wstar2, rstar * (1 - wstar2))
      }
      else if (marg == "Rd") {
        xstar <- pol2cart.L1Rd(r = rstar, w = wstar2)
      }
      return(list(xstar = xstar, iw = iw))
    }
  })
  return(sims)
}
