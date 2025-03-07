library(MCMCvis)
library(posterior)
library(cmdstanr)
library(dplyr)
library(tidyr)
library(scales)
source("egpd_functions.R")


fwi <- RcppSimdJson::fload("data/erc_fwi.json")$fwi
library(extRemes)
thres95 <- quantile(fwi, 0.95)
evd_fit <- fevd(x = fwi, threshold = thres95, type = "GP", time.units = "days")
t(as.matrix(evd_fit$results$par))

## Mixture of EGPD with exponential
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

csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_expo_negxi_take2_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)
# MCMCtrace(fitmcmc, ind = TRUE, open_pdf = TRUE)

fit_df <- fit |> as_draws_df() |>
  select(-c(lp__, .chain, .iteration,.draw)) |>
  select(-contains("_prime"))

egpd_params_post_expo <- fit_df |> select("kappa", "xi", "sigma_egpd") |> rename(sigma = sigma_egpd) |> colMeans()
expo_params_post_expo <- fit_df |> select("scale") |> colMeans()
mix_probs_post_expo <- fit_df |> select(contains("theta")) |> rename(egpd_prob = "theta[1]", norm_prob = "theta[2]") |> colMeans()

fwi_unif <- mix_cdf_expo(fwi, mix_probs_post_expo, egpd_params_post_expo, expo_params_post_expo)
hist(fwi_unif, freq=FALSE, breaks = 35)

## Mixture of EGPD with truncated exponential
dtexp <- function(x, xmax, rate) {
  if(x > xmax) {
    return(0)
  } else {
    numer <- dexp(x, rate)
    denom <- pexp(xmax, rate)
    return(numer/denom)
  }
}

# x <- seq(0,120,length.out = 500)
# pdf <- rep(NA, length(x))
# for(i in seq_along(x)) {
#   pdf[i] <- dtexp(x[i], 60, expo_params_post_texpo["rate"])
# }
# 
# plot(x, pdf)
# dtexp_v2 <- function(x, xmax, rate) {
#   numer <- dexp(x, rate)
#   denom <- pexp(xmax, rate)
#   return(numer/denom)
# }

ptexp <- function(data, xmax, rate) {
  if(data > xmax) {
    return(1)
  } else {
    numer <- pexp(data, rate)
    denom <- pexp(xmax, rate)
    return(numer/denom)
  }
}

mix_dens_trunc_expo <- function(data, trunc_pt, mix_probs, egpd_params, expo_params) {
  egpd_part <- g1_pdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
  expo_part <- rep(NA, length(data))
  for(i in seq_along(data)) {
    expo_part[i] <- dtexp(data[i], trunc_pt, expo_params["rate"])
  }
  return(mix_probs[1] * egpd_part + mix_probs[2] * expo_part)
}

# mix_dens_trunc_expo_v2 <- function(data, mix_probs, egpd_params, expo_params) {
#   egpd_part <- g1_pdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
#   expo_part <- dtexp_v2(data, 60, expo_params["rate"])
#   return(mix_probs[1] * egpd_part + mix_probs[2] * expo_part)
# }

mix_cdf_trunc_expo <- function(data, trunc_pt, mix_probs, egpd_params, expo_params) {
  egpd_part <- g1_cdf(data,sigma = egpd_params["sigma"], xi = egpd_params["xi"], kappa = egpd_params["kappa"])
  expo_part <- rep(NA, length(data))
  for(i in seq_along(data)) {
    expo_part[i] <- ptexp(data[i], trunc_pt, expo_params["rate"])
  }
  return(mix_probs[1] * egpd_part + mix_probs[2] * expo_part)
}

csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_trunc_expo_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)

# csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
#                        pattern = "fwi_mix_trunc_expo_higher_ub_\\d{1}.csv", full.names = TRUE)
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "fwi_mix_trunc_expo_\\d{1}.csv", full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)
# MCMCtrace(fitmcmc, ind = TRUE, open_pdf = TRUE)
# 
fit_df <- fit |> as_draws_df() |>
  select(-c(lp__, .chain, .iteration,.draw)) |>
  select(-contains("_prime"))

egpd_params_post_texpo60 <- fit_df |> select("kappa", "xi", "sigma_egpd") |> rename(sigma = sigma_egpd) |> colMeans()
expo_params_post_texpo60 <- fit_df |> select("rate") |> colMeans()
mix_probs_post_texpo60 <- fit_df |> select(contains("theta")) |> rename(egpd_prob = "theta[1]", norm_prob = "theta[2]") |> colMeans()

fwi_unif <- mix_cdf_trunc_expo(fwi, 80, mix_probs_post_texpo, egpd_params_post_texpo, expo_params_post_texpo)
hist(fwi_unif, freq = FALSE, breaks = 25)

hist(fwi, freq = FALSE, breaks = 45)
curve(mix_dens_trunc_expo(x, 60, mix_probs_post_texpo60, egpd_params_post_texpo60, expo_params_post_texpo60), add = TRUE, col = "red", lwd = 3)
curve(mix_dens_trunc_expo(x, 80, mix_probs_post_texpo80, egpd_params_post_texpo80, expo_params_post_texpo80), add = TRUE, col = "purple", lwd = 3)
curve(mix_dens_trunc_expo(x, 100, mix_probs_post_texpo100, egpd_params_post_texpo100, expo_params_post_texpo100), add = TRUE, col = "blue", lwd = 3)
curve(mix_dens_trunc_expo(x, 120, mix_probs_post_texpo120, egpd_params_post_texpo120, expo_params_post_texpo120), add = TRUE, col = "orange", lwd = 3)
# curve(mix_dens_trunc_expo_v2(x, mix_probs_post_texpo, egpd_params_post_texpo, expo_params_post_texpo), add = TRUE, col = alpha("blue", 0.75), lwd = 3)
# curve(mix_dens_expo(x, mix_probs_post_expo, egpd_params_post_expo, expo_params_post_expo), add = TRUE, col = alpha("orange", 0.75), lwd = 3)


## exponential margins
fwi_mix60_expo <- qexp(mix_cdf_trunc_expo(fwi, 60, mix_probs_post_texpo60, egpd_params_post_texpo60, expo_params_post_texpo60))
fwi_mix80_expo <- qexp(mix_cdf_trunc_expo(fwi, 80, mix_probs_post_texpo80, egpd_params_post_texpo80, expo_params_post_texpo80))
fwi_mix100_expo <- qexp(mix_cdf_trunc_expo(fwi, 100, mix_probs_post_texpo100, egpd_params_post_texpo100, expo_params_post_texpo100))
fwi_mix120_expo <- qexp(mix_cdf_trunc_expo(fwi, 120, mix_probs_post_texpo120, egpd_params_post_texpo120, expo_params_post_texpo120))

# qqplot(fwi_expo_mix, rexp(3000))
qqplot(fwi_texpo_mix, rexp(2945))
# qqplot(fwi_texpo_ub_mix, rexp(3000))

csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = "redstone_g1_\\d{1}.csv", full.names = TRUE)
fit_df <- as_cmdstan_fit(csvfiles)$draws() |> as_draws_df() |> select("xi[1]", "kappa[1]", "sigma[1]") |>
  rename(kappa = "kappa[1]", xi = "xi[1]", sigma = "sigma[1]") |>
  apply(MARGIN = 2, median)

erc <- RcppSimdJson::fload("data/erc_fwi.json")$erc
hist(erc, freq = FALSE)
erc_unif <- g1_cdf(erc, sigma = fit_df["sigma"], xi = fit_df["xi"], kappa = fit_df["kappa"])
hist(erc_unif, freq = FALSE, breaks =25)
curve(g1_pdf(x, sigma = fit_df["sigma"], xi = fit_df["xi"], kappa = fit_df["kappa"]), add = TRUE)

erc_expo <- qexp(g1_cdf(erc, sigma = fit_df["sigma"], xi = fit_df["xi"], kappa = fit_df["kappa"]))
qqplot(erc_expo, rexp(length(erc_expo)))

n <- length(erc)

plot(erc_expo/log(n), fwi_mix60_expo/log(n), pch = 20, col = alpha("purple", 0.5))
points(erc_expo/log(n), fwi_mix120_expo/log(n), pch = 20, col = alpha("blue", 0.5))
points(erc_expo/log(n), fwi_mix100_expo/log(n), pch = 20, col = alpha("red", 0.5))
points(erc_expo/log(n), fwi_mix80_expo/log(n), pch = 20, col = alpha("orange", 0.5))


exp_ERC_rank <- qexp(rank(erc)/(length(erc) + 1))
exp_FWI_rank <- qexp(rank(fwi)/(length(fwi) + 1))
points(exp_ERC_rank/log(2945), exp_FWI_rank/log(2945), pch = 20, col = alpha("purple", 0.5))

# proceeding with UB = 100 for FWI
indices_expo <- cbind(x = erc_expo, y = fwi_mix100_expo) |> as_tibble() |>
  mutate(r = x + y,
         w1 = x / r,
         w2 = y / r)
grab_top_n <- function(cloud_tib, n0 = 1, N = 5000) {
  tau <- (N-n0)/N
  q1 <- quantile(cloud_tib$x, tau)
  q2 <- quantile(cloud_tib$y, tau)
  q <- max(q1, q2)
  idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  eps <- 0.001
  while (length(idx) > n0) {
    q <- q + eps
    idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  }
  cloud_tib <- cloud_tib |> 
    mutate(r0_w = ifelse(w1 > 0.5, q/w1, q/w2),
           x_lb = ifelse(w1 < 0.5, q, q*y/x),
           y_lb = ifelse(w1 > 0.5, q, q*x/y),
           high = as.factor(case_when(y >= q | x >= q ~ 1,
                                      .default = 0)))
  return(list(q = q,
              idx = idx,
              n0 = length(idx),
              N = N,
              R = cloud_tib$r,
              W = cloud_tib$w1,
              W2 = cloud_tib$w2,
              r0_w = cloud_tib$r0_w,
              cloud_tib = cloud_tib))
  
}

fire_data_list <- grab_top_n(indices_expo, n0 = 150, N = nrow(indices_expo))
plot(indices_expo$x, indices_expo$y, col = alpha("blue", 0.5), pch = 20)
points(indices_expo$x[fire_data_list$idx], indices_expo$y[fire_data_list$idx], pch = 20, col = alpha("red", 0.5))

cols <- c("lightblue", "blue", "red")
fire_data_list$cloud_tib |>
  ggplot(aes(x/log(n), y/log(n), color = high)) + 
  geom_point(alpha=0.8) +
  geom_line(aes(x = y_lb/log(n), y = x_lb/log(n), color = 'red'), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  # scale_x_continuous(expand = c(0,0.01)) +
  # scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = expansion(mult = c(0,0.025))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.025))) +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  xlab(expression("ERC"/"log(n)")) + ylab(expression("FWI"/"log(n)"))

fire_data_list$cloud_tib |>
  ggplot(aes(w1, r/log(n), color = high)) + 
  geom_point(alpha=0.8) +
  geom_line(aes(x = x_lb/(x_lb + y_lb), y = (x_lb + y_lb)/log(n), color = 'red'), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.position = "none") +
  xlab(expression("W")) + ylab(expression("R"/"log(n)"))

qs::qsave(fire_data_list, file = "data/fire_data_expo.qs")
