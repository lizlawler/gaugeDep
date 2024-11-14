library(nimble)
library(tidyverse)

gauss_gauge <- nimbleFunction(
  run = function(W = double(0), dep = double(0)) {
    w2 <- 1 - W
    result <- (W + w2 - 2 * dep * sqrt(W * w2)) / (1 - dep^2)
    returnType(double(0))
    return(result)
  }
)

dangular <- nimbleFunction(
  run = function(x = double(0), pars = double(0),
                 vol_est = double(0),
                 dim = integer(0, default = 2),
                 log = integer(0, default = 1)) {
    returnType(double(0))
    logProb <- -dim * log(gauss_gauge(x, pars)) - (log(dim) + log(vol_est))
    if(log) return(logProb)
    else return(exp(logProb))
  })

est_vol <- nimbleFunction(
  run = function(sum_term = double(1), sqrt_term = double(1), pars = double(0)) {
    gx <- (sum_term - 2 * pars * sqrt_term) / (1 - pars^2)
    vol_est <- mean(gx <= 1)
    returnType(double(0))
    return(vol_est)
  }
)

ang_vol_code <- nimbleCode({
  # priors
  dep_w ~ dunif(0, 1)
  
  L_volume <- est_vol(sum_term[1:K], sqrt_term[1:K], dep_w)
  for(i in 1:N) {
    w[i] ~ dangular(pars = dep_w, vol_est = L_volume, dim = 2)
  }
})

grid <- expand.grid(seq(0, 1, length.out=100), seq(0, 1, length.out=100))
sum_grid <- grid[,1] + grid[,2]
sqrt_grid <- sqrt(grid[,1] * grid[,2])
data <- RcppSimdJson::fload("data/gauss/high_1.json")
angle_data <- list(
  w = data$W, 
  sum_term = sum_grid, 
  sqrt_term = sqrt_grid
)
angle_constants <- list(
  N = length(data$W),
  K = length(sum_grid)
)
angle_inits <- list(dep_w = runif(1))
angle_model <- nimbleModel(ang_vol_code, angle_constants, angle_data, angle_inits)
cmodel <- compileNimble(angle_model, resetFunctions = TRUE)
conf_model <- configureMCMC(angle_model)
modelMCMC <- buildMCMC(conf_model)

cmodelMCMC <- compileNimble(modelMCMC, project = angle_model)
model_run <- runMCMC(cmodelMCMC, niter = 20000, summary = TRUE, nburnin = 4000, nchains = 2)

