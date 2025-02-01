args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
gauge <- args[3]
likelihood <- args[4]

library(qs)
library(dplyr)
library(tidyr)

extract_median_params_radial <- function(dep_type, dep_level, gauge, data_num, likelihood) {
  params <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs",
                          dep_type, gauge, likelihood, dep_level, data_num))$samples |> as_tibble()
  return(params |> 
           select(-sigma_m) |>
           apply(MARGIN = 2, FUN = median) |> 
           t() |>
           as_tibble() |>
           mutate(dataset = data_num))
}


create_tib_med_params_radial <- function(dep_type, dep_level, gauge, likelihood) {
  data_num <- 1:100
  tib_med_params <- sapply(data_num, 
                           function(x) extract_median_params_radial(dep_type = dep_type, 
                                                                    dep_level = dep_level, 
                                                                    gauge = gauge, 
                                                                    data_num = x, 
                                                                    likelihood = likelihood), 
                           simplify = FALSE)
  return(tib_med_params |> bind_rows())
}

med_params <- create_tib_med_params_radial(dep_type = dep_type, 
                                           dep_level = dep_level, 
                                           gauge = gauge, 
                                           likelihood = likelihood)

filepath <- sprintf("fits_and_weights/med_params_joint/%s_%s_%s_radial.qs", 
                    dep_type, dep_level, likelihood)
qsave(med_params, filepath)
print("Posterior medians of radial parameters have been successfully saved")
