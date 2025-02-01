library(cmdstanr)
library(tidyverse)
library(MCMCvis)
library(posterior)
# library(purrr)
source("egpd_functions.R")

## read in stan results
read_fit <- function(location) {
  csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                         pattern = paste0(location, "_g1_exp_\\d{1}.csv"),
                         full.names = TRUE)
  return(read_cmdstan_csv(csvfiles, variables = c("xi", "kappa", "sigma"))$post_warmup_draws)
}

extract_med_params <- function(cmdstan_obj) {
  temp <- cmdstan_obj |> as_draws_df() |> 
    rename(draw = ".draw") |>
    select(!contains(c("log", "lp", "chain", "iter"))) |>
    pivot_longer(cols = -"draw") |> 
    separate_wider_delim(name, delim = "\\[", names = c("param", "index")) |>
    mutate(index = as.numeric(gsub("]\\", "", index)),
           index = case_when(index == 1 ~ "erc",
                             index == 2 ~ "fwi")) |> 
    group_by(param, index) |>
    summarize(post_med = median(value)) |> ungroup() |>
    pivot_wider(names_from = param, values_from = post_med)
}

exp_fit <- read_fit("redstone") |> as_draws_df() |> 
  rename(draw = ".draw") |>
  select(!contains(c("log", "lp", "chain", "iter"))) |>
  pivot_longer(cols = -"draw") |> 
  separate_wider_delim(name, delim = "[", names = c("param", "index")) |>
  mutate(index = as.numeric(gsub("]", "", index)),
         index = case_when(index == 1 ~ "erc",
                           index == 2 ~ "fwi")) |> 
  group_by(param, index) |>
  summarize(post_med = median(value)) |> ungroup() |>
  pivot_wider(names_from = param, values_from = post_med)

csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                                      pattern = paste0("redstone", "_g1_exp_\\d{1}.csv"),
                                      full.names = TRUE)
exp_fit_mcmc <- as_cmdstan_fit(csvfiles)

mcmc <- coda::as.mcmc.list(test_fit)

fitmcmc <- as_mcmc.list(exp_fit_mcmc)

## uncomment the below if directly interfacing with RStudio server on a remote machine
# options(bitmapType='cairo')

print("Creating traceplot of rhos...")
MCMCtrace(fitmcmc,
          params = c('xi', 'kappa', 'sigma'), 
          ind = TRUE, 
          open_pdf = TRUE)

test_params <- extract_med_params(test_fit)
og_data <- RcppSimdJson::fload("data/erc_fwi_og.json")

marg_transform <- function(data, egpd_family, ...) {
  egpd_cdf <- get(paste0(egpd_family, "_cdf"))
  exp_data <- qexp(egpd_cdf(data, ...))
  return(exp_data)
}

g1_marg_transform <- function(data, med_params_tib) {
  return(marg_transform(data, "g1", 
                        kappa = med_params_tib$kappa,
                        sigma = med_params_tib$sigma,
                        xi = med_params_tib$xi))
}

g2_marg_transform <- function(data, med_params_tib) {
  return(marg_transform(data, "g2", 
                        kappa1 = med_params_tib$kappa1,
                        kappa2 = med_params_tib$kappa2,
                        p = med_params_tib$probs,
                        sigma = med_params_tib$sigma,
                        xi = med_params_tib$xi))
}

g3_marg_transform <- function(data, med_params_tib) {
  return(marg_transform(data, "g3", 
                        sigma = med_params_tib$sigma,
                        xi = med_params_tib$xi,
                        delta = med_params_tib$delta))
}

g4_marg_transform <- function(data, med_params_tib) {
  return(marg_transform(data, "g4", 
                        kappa = med_params_tib$kappa,
                        sigma = med_params_tib$sigma,
                        xi = med_params_tib$xi,
                        delta = med_params_tib$delta))
}

read_extract_transform <- function(location, dataset, egpd_family) {
  fit <- read_fit(location, dataset, egpd_family)
  post_params <- extract_med_params(fit)
  fwi_params <- post_params |> filter(index == "fwi")
  erc_params <- post_params |> filter(index == "erc")
  data <- RcppSimdJson::fload(paste0("data/erc_fwi_", dataset, ".json"))
  egpd_transform <- get(paste0(egpd_family, "_marg_transform"))
  
  fwi_exp <- egpd_transform(data$fwi, fwi_params)
  erc_exp <- egpd_transform(data$erc, erc_params)
  return(tibble(erc = erc_exp, fwi = fwi_exp))
}

erc_unif <- g1_cdf(erc, test_fit$sigma[1], test_fit$xi[1], test_fit$kappa[1])
hist(erc_unif)

fwi_unif <- g1_cdf(fwi, test_fit$sigma[2], test_fit$xi[2], test_fit$kappa[2])
hist(fwi_unif)



points(erc, erc_exp)
test_g1 <- read_extract_transform("redstone", "og", "g1")
plot(test_g1/log(nrow(test_g1)), pch = 20, xlim = c(0,1), ylim = c(0,1))

test_g2 <- read_extract_transform("redstone", "og", "g2")
points(test_g2/log(nrow(test_g2)), pch = 20, col = "blue")

test_g3 <- read_extract_transform("redstone", "og", "g3")
points(test_g3/log(nrow(test_g3)), pch = 20, col = "red")

test_g4 <- read_extract_transform("redstone", "og", "g4")
points(test_g4/log(nrow(test_g4)), pch = 20, col = "green")

fwi_rank <- qexp(rank(og_data$fwi)/(length(og_data$fwi) + 1))
erc_rank <- qexp(rank(og_data$erc)/(length(og_data$erc) + 1))
fwi_erc_rank <- tibble(erc = erc_rank, fwi = fwi_rank)

plot(fwi_erc_rank/log(nrow(fwi_erc_rank)), pch = 20, col = "orange", xlim = c(0,1), ylim = c(0,1))

exp_ERC_rank <- qexp(rank(erc_data)/(length(erc_data) + 1))
exp_FWI_rank <- qexp(rank(fwi_data)/(length(fwi_data) + 1))

fwi_exp <- marg_transform(og_data$fwi, "g1", unlist(c(test_params[2,2:4])))
erc_exp <- g1_marg_transform(og_data$erc, test_params[1,])

test_params <- extract_med_params(test)
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/", pattern = ".csv",full.names = TRUE)
fit <- as_cmdstan_fit(csvfiles)
fit_mcmc <- as_mcmc.list(fit)

MCMCtrace(fit_mcmc, ind = TRUE)
extract_med_params <- function()
post_params <- fit |> as_draws_df() |> 
  rename(draw = ".draw") |>
  select(!contains(c("log", "lp", "chain", "iter"))) |>
  pivot_longer(cols = -"draw") |> 
  separate_wider_delim(name, delim = "[", names = c("param", "index")) |>
  mutate(index = as.numeric(gsub("]", "", index)),
         index = case_when(index == 1 ~ "erc",
                           index == 2 ~ "fwi")) |> 
  group_by(param, index) |>
  summarize(post_med = median(value)) |> ungroup()

erc_data <- redstone_summer$ERC[-idx]
fwi_data <- redstone_summer$FWI[-idx]
erc_params <- poster_params |> filter(index == "erc") |> select(-index) |> pivot_wider(names_from = param, values_from = post_med)
fwi_params <- poster_params |> filter(index == "fwi") |> select(-index) |> pivot_wider(names_from = param, values_from = post_med)
expERC <- qexp(egpd_cdf(erc_data, sigma = erc_params$sigma, xi = erc_params$xi, kappa = erc_params$kappa))
expFWI <- qexp(egpd_cdf(fwi_data, sigma = fwi_params$sigma, xi = fwi_params$xi, kappa = fwi_params$kappa))
plot(expERC/log(length(erc_data)), expFWI/log(length(erc_data)), pch = 20, xlim = c(0,1), ylim = c(0,1))
points(exp_ERC_rank/log(length(erc_data)), exp_FWI_rank/log(length(erc_data)), pch = 20, col = "blue")
qqplot(expERC, exp_ERC_rank)
qqplot(expFWI, exp_FWI_rank)
abline(a = 0, b = 1)
erc_dens <- egpd_pdf(erc_data, sigma = erc_params$sigma, xi = erc_params$xi, kappa = erc_params$kappa)
fwi_dens <- egpd_pdf(fwi_data, sigma = fwi_params$sigma, xi = fwi_params$xi, kappa = fwi_params$kappa)
erc_full <- tibble(vals = erc_data, density = erc_dens)
erc_full |> ggplot(aes(x = vals)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 boundary = 0, color = "black", fill = "lightgrey") +
  geom_line(aes(x = vals, y = density), linewidth = 1, alpha = 0.8)