library(nimble)
library(tidyverse)

gauss_gauge <- nimbleFunction(
  run = function(W = double(1), dep = double(0)) {
    w2 <- 1 - W
    result <- (W + w2 - 2 * dep * sqrt(W * w2)) / (1 - dep^2)
    returnType(double(1))
    return(result)
  }
)

est_vol <- nimbleFunction(
  run = function(sum_term = double(1), sqrt_term = double(1), pars = double(0)) {
    gx <- (sum_term - 2 * pars * sqrt_term) / (1 - pars^2)
    vol_est <- mean(gx <= 1)
    returnType(double(0))
    return(vol_est)
  }
)

angular_loglik <- nimbleFunction(
  run = function(W = double(1), pars = double(0),
                 sum_term = double(1), sqrt_term = double(1), dim = integer(0)) {
    N <- length(W)
    L_volume <- est_vol(sum_term, sqrt_term, pars)
    loglik <- -(dim * sum(log(gauss_gauge(W, pars)))) - (N * (log(dim) + log(L_volume)))
    returnType(double(0))
    return(loglik)
  }
)

# dangular <- nimbleFunction(
#   run = function(x = double(0), pars = double(0),
#                  vol_est = double(0),
#                  dim = integer(0, default = 2),
#                  log = integer(0, default = 1)) {
#     returnType(double(0))
#     logProb <- -dim * log(gauss_gauge(x, pars)) - (log(dim) + log(vol_est))
#     if(log) return(logProb)
#     else return(exp(logProb))
#   })
# 
# dradii <- nimbleFunction(
#   run = function(x = double(0), threshold = double(0), 
#                  par1 = double(0), par2 = double(0),
#                  log = integer(0, default = 1)) {
#     returnType(double(0))
#     if(x < threshold) {
#       logProb <- pgamma(threshold, shape = par1, rate = par2, log.p = TRUE)
#       if(log) return(log)
#       else return(exp(log))
#     } else {
#       logProb <- dgamma(x, shape = par1, rate = par2, log = TRUE)
#       if(log) return(log)
#       else return(exp(log))
#     }
#   }
# )

radial_cens_loglik <- nimbleFunction(
  run = function(radii = double(1), pars = double(1),
                 threshold = double(1), W = double(1)) {
    N <- length(radii)
    alpha <- pars[1]
    dep <- pars[2]
    beta <- gauss_gauge(W, dep)
    loglik <- 0.0

    for(i in 1:N) {
      if(radii[i] < threshold[i]) {
        loglik <- loglik + pgamma(threshold[i], alpha, rate = beta[i], log = TRUE)
      } else {
        loglik <- loglik + dgamma(radii[i], alpha, rate = beta[i], log = TRUE)
      }
    }

    returnType(double(0))
    return(loglik)
  }
)

joint_loglik <- nimbleFunction(
  run = function(R = double(1), W = double(1), params = double(1), threshold = double(1),
                 sum_term = double(1), sqrt_term = double(1), dim = integer(0)) {

    radial_ll <- radial_cens_loglik(R, c(params[1], params[2]), threshold, W)
    angular_ll <- angular_loglik(W, params[3], sum_term, sqrt_term, dim)
    total_ll <- radial_ll + angular_ll

    returnType(double(0))
    return(total_ll)
  }
)

joint_code <- nimbleCode({
  # priors
  alpha ~ dgamma(4, 2)
  dep_r ~ dunif(0, 1)
  dep_w ~ dunif(0, 1)
  
  # # estimate volume for angular density
  # L_volume <- est_vol(sum_term[1:K], sqrt_term[1:K], dep_w)
  # 
  # # likelihood
  # for(i in 1:N) {
  #   w[i] ~ dangular(dep_w, L_volume)
  #   r[i] ~ dradii(threshold[i], alpha, dep_r)
  # }
  loglik <- joint_loglik(R[1:N], W[1:N], c(alpha, dep_r, dep_w),
                         threshold[1:N], sum_term[1:K], sqrt_term[1:K], 2)
})

grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))
sum_grid <- grid[,1] + grid[,2]
sqrt_grid <- sqrt(grid[,1] * grid[,2])
data <- RcppSimdJson::fload("data/gauss/high_1.json")
joint_data <- list(
  R = data$R, 
  W = data$W, 
  threshold = data$r0_w,
  sum_term = sum_grid, 
  sqrt_term = sqrt_grid
)
joint_constants <- list(
  N = length(data$R),
  K = length(sum_grid)
)
joint_inits <- list(alpha = rgamma(1, 4, 2), dep_r = runif(1), dep_w = runif(1))
joint_model <- nimbleModel(joint_code, joint_constants, joint_data, joint_inits)
cmodel <- compileNimble(joint_model, resetFunctions = TRUE)
conf_model <- configureMCMC(joint_model)
conf_model$removeSamplers()
conf_model$addSampler(target = c('alpha', 'dep_r', 'dep_w'), 'AF_slice')
modelMCMC <- buildMCMC(conf_model)

cmodelMCMC <- compileNimble(modelMCMC, project = joint_model, resetFunctions = TRUE)
model_run <- runMCMC(cmodelMCMC, niter = 20000, summary = TRUE, nburnin = 4000, nchains = 2)
