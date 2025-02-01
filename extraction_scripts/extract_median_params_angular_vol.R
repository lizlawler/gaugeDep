args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
gauge <- args[3]

library(qs)
library(dplyr)
library(tidyr)

extract_median_params_ang <- function(dep_type, dep_level, gauge, data_num) {
  params <- qread(sprintf("samplers/rcpp/angular_vol_mcmc_fits/%s/%s_%s_%s.qs",
                          dep_type, gauge, dep_level, data_num))$samples |> as_tibble()
  return(params |> 
           select(-sigma_m) |>
           apply(MARGIN = 2, FUN = median) |> 
           t() |>
           as_tibble() |>
           mutate(dataset = data_num))
}


create_tib_med_params_ang <- function(dep_type, dep_level, gauge) {
  data_num <- 1:100
  tib_med_params <- sapply(data_num, 
                           function(x) extract_median_params_ang(dep_type = dep_type, 
                                                                 dep_level = dep_level, 
                                                                 gauge = gauge, 
                                                                 data_num = x), 
                           simplify = FALSE)
  return(tib_med_params |> bind_rows())
}

med_params <- create_tib_med_params_ang(dep_type = dep_type, 
                                        dep_level = dep_level, 
                                        gauge = gauge)

filepath <- sprintf("fits_and_weights/med_params_joint/%s_%s_%s_ang_vol.qs", 
                    dep_type, dep_level, gauge)
qsave(med_params, filepath)
print("Posterior medians of star angular parameters have been successfully saved")
