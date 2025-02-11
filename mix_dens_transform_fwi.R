library(MCMCvis)
library(posterior)
library(cmdstanr)
library(dplyr)
library(tidyr)
library(scales)
source("egpd_functions.R")

## Mixture of EGPD with normal, xi > 0
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)
# MCMCtrace(fitmcmc, ind = TRUE, open_pdf = TRUE)

fit_df <- fit |> as_draws_df() |>
  select(-c(lp__, .chain, .iteration,.draw)) |>
  select(-contains("_prime"))

egpd_params_post <- fit_df |> select("kappa", "xi", "sigma_egpd") |> rename(sigma = sigma_egpd) |> colMeans()
normal_params_post <- fit_df |> select("sigma_norm", "mu") |> rename(sigma = sigma_norm) |> colMeans()
mix_probs_post <- fit_df |> select(contains("theta")) |> rename(egpd_prob = "theta[1]", norm_prob = "theta[2]") |> colMeans()

mix_dens <- function(data, mix_probs, egpd_params, normal_params) {
  egpd_part <- g1_pdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
  normal_part <- dnorm(data, mean = normal_params["mu"], sd = normal_params["sigma"])
  return(mix_probs[1] * egpd_part + mix_probs[2] * normal_part)
}

mix_cdf <- function(data, mix_probs, egpd_params, normal_params) {
  egpd_part <- g1_cdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
  normal_part <- pnorm(data, mean = normal_params["mu"], sd = normal_params["sigma"])
  return(mix_probs[1] * egpd_part + mix_probs[2] * normal_part)
}

fwi <- RcppSimdJson::fload("data/erc_fwi.json")$fwi
fwi_unif <- mix_cdf(fwi, mix_probs_post, egpd_params_post, normal_params_post)
hist(fwi, freq = FALSE, breaks = 40)
curve(mix_dens(x, mix_probs_post, egpd_params_post, normal_params_post), add = TRUE)

## Mixture of EGPD with exponential, xi > 0
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_expo_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)
# MCMCtrace(fitmcmc, ind = TRUE, open_pdf = TRUE)

fit_df <- fit |> as_draws_df() |>
  select(-c(lp__, .chain, .iteration,.draw)) |>
  select(-contains("_prime"))

egpd_params_post_expo <- fit_df |> select("kappa", "xi", "sigma_egpd") |> rename(sigma = sigma_egpd) |> colMeans()
expo_param_post <- fit_df |> select("scale") |> colMeans()
mix_probs_post_expo <- fit_df |> select(contains("theta")) |> rename(egpd_prob = "theta[1]", norm_prob = "theta[2]") |> colMeans()

mix_dens_expo <- function(data, mix_probs, egpd_params, expo_params) {
  egpd_part <- g1_pdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
  expo_part <- dexp(data, expo_params["scale"])
  return(mix_probs[1] * egpd_part + mix_probs[2] * expo_part)
}

mix_cdf_expo <- function(data, mix_probs, egpd_params, expo_params) {
  egpd_part <- g1_cdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
  expo_part <- pexp(data, expo_params["scale"])
  return(mix_probs[1] * egpd_part + mix_probs[2] * expo_part)
}

fwi_unif <- mix_cdf_expo(fwi, mix_probs_post_expo, egpd_params_post_expo, expo_param_post)
hist(fwi, freq = FALSE, breaks = 50)
curve(mix_dens(x, mix_probs_post, egpd_params_post, normal_params_post), add = TRUE, col = "blue")
curve(mix_dens_expo(x, mix_probs_post_expo, egpd_params_post_expo, expo_param_post), add = TRUE, col = "red")

library(extRemes)
thres95 <- quantile(fwi, 0.95)
evd_fit <- fevd(x = fwi, threshold = thres95, type = "GP", time.units = "days")
t(as.matrix(evd_fit$results$par))

## Mixture of EGPD with exponential, xi < 0
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_expo_negxi_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)
MCMCtrace(fitmcmc, ind = TRUE, open_pdf = TRUE)

fit_df <- fit |> as_draws_df() |>
  select(-c(lp__, .chain, .iteration,.draw)) |>
  select(-contains("_prime"))

egpd_params_post_expo_negxi <- fit_df |> select("kappa", "xi", "sigma_egpd") |> rename(sigma = sigma_egpd) |> colMeans()
expo_params_post_negxi <- fit_df |> select("scale") |> colMeans()
mix_probs_post_negxi <- fit_df |> select(contains("theta")) |> rename(egpd_prob = "theta[1]", norm_prob = "theta[2]") |> colMeans()

fwi_unif <- mix_cdf_expo(fwi, mix_probs_post_negxi, egpd_params_post_expo_negxi, expo_params_post_negxi)
fwi <- RcppSimdJson::fload("data/erc_fwi.json")$fwi
hist(fwi, freq = FALSE, breaks = 45)
curve(mix_dens(x, mix_probs_post, egpd_params_post, normal_params_post), add = TRUE, col = "red", lwd = 3)
curve(mix_dens_expo(x, mix_probs_post_expo, egpd_params_post_expo, expo_param_post), add = TRUE, col = "blue", lwd = 3)
curve(mix_dens_expo(x, mix_probs_post_negxi, egpd_params_post_expo_negxi, expo_params_post_negxi), add = TRUE, col = "orange", lwd = 3)


## exponential margins
fwi_expo_norm_mix <- qexp(mix_cdf(fwi, mix_probs_post, egpd_params_post, normal_params_post))
fwi_expo_expo_mix <- qexp(mix_cdf_expo(fwi, mix_probs_post_expo, egpd_params_post_expo, expo_param_post))
fwi_expo_expo_negxi_mix <- qexp(mix_cdf_expo(fwi, mix_probs_post_negxi, egpd_params_post_expo_negxi, expo_params_post_negxi))

csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "redstone_g1_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)

fit_df <- fit$draws() |> as_draws_df() |> select("xi[1]", "kappa[1]", "sigma[1]") |>
  rename(kappa = "kappa[1]", xi = "xi[1]", sigma = "sigma[1]") |>
  apply(MARGIN = 2, median)

erc <- RcppSimdJson::fload("data/erc_fwi.json")$erc
hist(erc, freq = FALSE)
curve(g1_pdf(x, sigma = fit_df["sigma"], xi = fit_df["xi"], kappa = fit_df["kappa"]), add = TRUE)

erc_expo <- qexp(g1_cdf(erc, sigma = fit_df["sigma"], xi = fit_df["xi"], kappa = fit_df["kappa"]))
plot(erc_expo/log(2945), fwi_expo_norm_mix/log(2945), pch = 20, col = alpha("blue", 0.5), ylim = c(0,1))
points(erc_expo/log(2945), fwi_expo_expo_mix/log(2945), pch = 20, col = alpha("red", 0.5))
plot(erc_expo, fwi_expo_expo_negxi_mix, pch = 20, col = alpha("orange", 0.5), xlim = c(0, 14), ylim = c(0,14))


exp_ERC_rank <- qexp(rank(erc)/(length(erc) + 1))
exp_FWI_rank <- qexp(rank(fwi)/(length(fwi) + 1))
points(exp_ERC_rank/log(2945), exp_FWI_rank/log(2945), pch = 20, col = alpha("purple", 0.5))

## Mixture of EGPD with half-norm, xi < 0
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_halfnorm_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)
MCMCtrace(fitmcmc, ind = TRUE, open_pdf = TRUE)

fit_df <- fit |> as_draws_df() |>
  select(-c(lp__, .chain, .iteration,.draw)) |>
  select(-contains("_prime"))

egpd_params_norm_negxi <- fit_df |> select("kappa", "xi", "sigma_egpd") |> rename(sigma = sigma_egpd) |> colMeans()
norm_params_negxi <- fit_df |> select("sigma_norm", "mu") |> rename(sigma = sigma_norm) |> colMeans()
mix_probs_norm_negxi <- fit_df |> select(contains("theta")) |> rename(egpd_prob = "theta[1]", norm_prob = "theta[2]") |> colMeans()

half_norm_pdf <- function(x, sigma) {
  lpdf <- 0.5 * log(2) - log(sigma) - log(pi) - x^2 / (2 * sigma^2)
  return(exp(lpdf))
}
# mix_dens_halfnorm <- function(data, mix_probs, egpd_params, half_norm_params) {
#   egpd_part <- g1_pdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
#   hn_part <- half_norm_pdf(data, half_norm_params["sigma"])
#   return(mix_probs[1] * egpd_part + mix_probs[2] * hn_part)
# }

hist(fwi,freq = FALSE, breaks = 45)
curve(half_norm_pdf(x, mix_probs_norm_negxi, egpd_params_norm_negxi, norm_params_negxi), add = TRUE, col = "blue", lwd = 3)

## Mixture of EGPD with exponential, xi < 0, take 2
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_expo_negxi_take2_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)
MCMCtrace(fitmcmc, ind = TRUE, open_pdf = TRUE)

fit_df <- fit |> as_draws_df() |>
  select(-c(lp__, .chain, .iteration,.draw)) |>
  select(-contains("_prime"))

egpd_params_post_expo_negxi <- fit_df |> select("kappa", "xi", "sigma_egpd") |> rename(sigma = sigma_egpd) |> colMeans()
expo_params_post_negxi <- fit_df |> select("scale") |> colMeans()
mix_probs_post_negxi <- fit_df |> select(contains("theta")) |> rename(egpd_prob = "theta[1]", norm_prob = "theta[2]") |> colMeans()

fwi_unif <- mix_cdf_expo(fwi, mix_probs_post_negxi, egpd_params_post_expo_negxi, expo_params_post_negxi)
hist(fwi, freq = FALSE, breaks = 45)
curve(mix_dens(x, mix_probs_post, egpd_params_post, normal_params_post), add = TRUE, col = "red", lwd = 3)
curve(mix_dens_expo(x, mix_probs_post_expo, egpd_params_post_expo, expo_param_post), add = TRUE, col = "blue", lwd = 3)
curve(mix_dens_expo(x, mix_probs_post_negxi, egpd_params_post_expo_negxi, expo_params_post_negxi), add = TRUE, col = "green", lwd = 3)
curve(mix_dens(x, mix_probs_norm_negxi, egpd_params_norm_negxi, norm_params_negxi), add = TRUE, col = "orange", lwd = 3)

