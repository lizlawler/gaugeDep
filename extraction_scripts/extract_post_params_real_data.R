library(qs)
library(dplyr)
library(tidyr)

# data_name <- "redstone"

extract_post_params_radial <- function(gauge, likelihood, data) {
  params <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/real_data/%s_%s_%s.qs",
                          data, gauge, likelihood))$samples |> as_tibble() |> 
           select(-sigma_m) |>
           apply(MARGIN = 2, FUN = median) |> 
           t() |>
           as_tibble()
  return(params)
}

extract_post_params_ang_star <- function(gauge, data) {
  params <- qread(sprintf("samplers/rcpp/ang_star_mcmc_fits/real_data/%s_%s.qs",
                          data, gauge))$samples |> as_tibble() |> 
           select(-sigma_m) |>
           apply(MARGIN = 2, FUN = median) |> 
           t() |>
           as_tibble()
  return(params)
}

extract_post_params_ang_mix <- function(data) {
  params <- qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/real_data/%s.qs",
                          data)) |>
    as_tibble() |>
    select(matches("probs|alphastar|betastar")) |>
           colMeans() |>
           t() |>
           as_tibble()
  return(params)
}

# extract_post_params_ang_mix(data_name)

