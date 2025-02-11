library(nimble)
library(tidyverse)
library(posterior)

sb_code <- nimbleCode({
  for(i in 1:(L-1)){
    v[i] ~ dbeta(1, alpha)
  }
  alpha ~ dgamma(2, 2)
  beta_tau ~ dexp(1/8)
  probs[1:L] <- stick_breaking(v[1:(L-1)])
  for(i in 1:L) {
    mustar[i] ~ dunif(0, 1)
    taustar[i] ~ dinvgamma(2, beta_tau)
    alphastar[i] <- mustar[i] * taustar[i]
    betastar[i] <- (1-mustar[i]) * taustar[i]
  }
  for(i in 1:N) {
    z[i] ~ dcat(probs[1:L])
    w[i] ~ dbeta(alphastar[z[i]], betastar[z[i]])
  }
})


sb_bern_code <- nimbleCode({
  for(i in 1:(L-1)){
    v[i] ~ dbeta(1, alpha)
  }
  alpha ~ dgamma(2, 2)
  probs[1:L] <- stick_breaking(v[1:(L-1)])
  for(j in 1:L) {
    alphastar[j] <- j
    betastar[j] <- L - j + 1
  }
  for(i in 1:N) {
    z[i] ~ dcat(probs[1:L])
    w[i] ~ dbeta(alphastar[z[i]], betastar[z[i]])
  }
})


data <- RcppSimdJson::fload("data/gauss/mid_5.json")
w <- data$W
idx <- data$idx
sb_data <- list(w = w[idx])
sb_constants <- list(N = length(sb_data$w), L=10)
sb_inits <- list(mustar = runif(sb_constants$L, 0, 1), 
                 taustar = rinvgamma(sb_constants$L, 1, 1), 
                 z = sample(1:7, size = sb_constants$N, replace = TRUE),
                 v  = rbeta(sb_constants$L, 1, 1),
                 alpha = 1,
                 beta_tau = 8)
sb_model <- nimbleModel(sb_code, sb_constants, sb_data, sb_inits)
cmodel <- compileNimble(sb_model)
conf_model <- configureMCMC(sb_model)
conf_model$addMonitors("alphastar", "betastar", "taustar", "probs", "z")
modelMCMC <- buildMCMC(conf_model)

cmodelMCMC <- compileNimble(modelMCMC, project = sb_model)
model_run <- runMCMC(cmodelMCMC, niter = 15000, nburnin = 5000, thin = 5, inits = sb_inits, nchains = 1)

mix_dens <- function(w, chain_of_params) {
  alphas <- chain_of_params[,grepl("alphastar", colnames(chain_of_params))]
  betas <- chain_of_params[,grepl("betastar", colnames(chain_of_params))]
  weights <- chain_of_params[,grepl("probs", colnames(chain_of_params))]
  
  n_clust <- ncol(weights)
  n_iter <- nrow(weights)
  n_angle <- length(w)
  dens <- matrix(0.0, nrow = n_iter, ncol = n_angle)
  
  for(j in 1:n_iter) {
    for(k in 1:n_clust) {
      dens[j, ] <- dens[j, ] + weights[j, k] * dbeta(w, alphas[j, k], betas[j, k])
    }
  }
  return(colMeans(dens))
}

hist(w[idx], freq = FALSE, xlim = c(0,1), breaks = 25)
curve(mix_dens(x, model_run), add = TRUE)


data <- RcppSimdJson::fload("data/gauss/high_1.json")
w <- data$W
idx <- data$idx
sb_bern_data <- list(w = w[idx])
sb_bern_constants <- list(N = length(sb_bern_data$w), L=12)
sb_bern_inits <- list(z = sample(1:floor(sb_bern_constants$L/2), size = sb_bern_constants$N, replace = TRUE),
                 v  = rbeta(sb_bern_constants$L, 1, 1),
                 alpha = 1)
dim_list <- list(alphastar = sb_bern_constants$L, betastar = sb_bern_constants$L)
sb_bern_model <- nimbleModel(sb_bern_code, sb_bern_constants, sb_bern_data, sb_bern_inits, dimensions = dim_list)
cmodel <- compileNimble(sb_bern_model)
conf_model <- configureMCMC(sb_bern_model)
conf_model$addMonitors("alphastar", "betastar", "probs", "z")
modelMCMC <- buildMCMC(conf_model)

cmodelMCMC <- compileNimble(modelMCMC, project = sb_bern_model)
model_run <- runMCMC(cmodelMCMC, niter = 15000, nburnin = 5000, thin = 5, inits = sb_bern_inits, nchains = 1)


z <- model_run[,grepl(pattern = "z\\[", colnames(model_run))]
table(z[2000,])
post_params <- colMeans(model_run)

mix_dens <- function(w, mean_params) {
  alphas <- as.numeric(mean_params[grepl("alphastar", names(mean_params))])
  betas <- as.numeric(mean_params[grepl("betastar", names(mean_params))])
  weights <- as.numeric(mean_params[grepl("probs", names(mean_params))])
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}
test <- mix_dens(w[idx], post_params)

hist(w[idx], freq = FALSE, xlim = c(0,1), breaks = 25)
curve(mix_dens(x, post_params), add = TRUE, col = "red")
curve(dbeta(x, 2, 6), add = TRUE)
points(sb_data$w, w_dens_sb)

files <- list.files(path = paste0("stan/radial_angular/csv_fits/", "gauss", "/"),
                    pattern = paste0("high", "_", 1, "_\\d{1}.csv"), full.names = TRUE)
post_params <- as_cmdstan_fit(files) |> as_draws_df() |>
  select(any_of(contains(c("weights", "alpha","beta", "dep", ".draw", ".chain")))) |>
  group_by(.chain) |>
  summarize(across(everything(), mean)) |>
  select(-.draw) |>
  filter(.chain == 3) |> 
  select(-.chain) |> t() |> as.data.frame() |>
  rownames_to_column() |>
  as_tibble()

stan_wts <- post_params |> filter(grepl("weights", rowname)) |> select(V1) |> unlist() |> as.numeric()
stan_alphas <- post_params |> filter(grepl("alpha_ang", rowname)) |> select(V1) |> unlist() |> as.numeric()
stan_betas <- post_params |> filter(grepl("beta", rowname)) |> select(V1) |> unlist() |> as.numeric()

w_dens_stan <- mix_dens(sb_data$w, stan_wts, stan_alphas, stan_betas)
points(sb_data$w, w_dens_stan, col = "green")
