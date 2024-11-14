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

sb_data <- list(w = RcppSimdJson::fload("data/gauss/low_1.json")$W)
sb_constants <- list(N = length(sb_data$w), L=10)
sb_inits <- list(mustar = runif(sb_constants$L, 0, 1), 
                 taustar = rinvgamma(sb_constants$L, 1, 1), 
                 z = sample(1:10, size = sb_constants$N, replace = TRUE),
                 v  = rbeta(sb_constants$L, 1, 1),
                 alpha = 1,
                 beta_tau = 8)
sb_model <- nimbleModel(sb_code, sb_constants, sb_data, sb_inits)
cmodel <- compileNimble(sb_model, dirName = "nimble/")
conf_model <- configureMCMC(sb_model)
conf_model$addMonitors("alphastar", "betastar", "taustar", "probs", "z")
modelMCMC <- buildMCMC(conf_model)

cmodelMCMC <- compileNimble(modelMCMC, project = sb_model, dirName = "nimble/")
model_run <- runMCMC(cmodelMCMC, niter = 15000, inits = sb_inits, nchains = 1)
# means_medians_post <- model_run$summary |> as.data.frame() |> rownames_to_column() |> select(rowname, Mean, Median)

z <- model_run[,grepl(pattern = "z\\[", colnames(model_run))]
table(z[15000,])
model_run <- gauss_low1_sb_params
model_probs <- model_run[,grepl(pattern = "probs\\[", colnames(model_run))]
model_alphas <- model_run[,grepl(pattern = "alphastar\\[", colnames(model_run))]
model_mus <- model_run[,grepl(pattern = "mustar\\[", colnames(model_run))]
model_betas <- model_run[,grepl(pattern = "betastar\\[", colnames(model_run))]

sb_probs <- colMeans(model_probs[7500:15000,])
sb_alphas <- colMeans(model_alphas[7500:15000,])
sb_betas <- colMeans(model_betas[7500:15000,])

sb_mus <- colMeans(model_mus[7500:15000,])
sb_taus <- sb_alphas / sb_mus

sb_mus
sb_taus
# cx2 <- model_run$samples$chain2
# cx2_z <- cx2[,grepl(pattern = "z\\[", colnames(cx2))]
# table(cx2_z[15000,])
# 
# cx3 <- model_run$samples$chain3
# cx3_z <- cx3[,grepl(pattern = "z\\[", colnames(cx3))]
# table(cx3_z[15000,])

plot(cx2[,"probs[8]"], type = "l")
saveRDS(model_run, file = "nimble/gauss_low1_sb_params.rds")


mix_dens <- function(w, post_wts, post_alpha, post_beta) {
  n <- length(post_wts)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + post_wts[i] * dbeta(w, post_alpha[i], post_beta[i])
  }
  return(dens)
}

w_dens_sb <- mix_dens(sb_data$w, sb_probs, sb_alphas, sb_betas)
hist(sb_data$w, freq = FALSE, breaks = 35)
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
