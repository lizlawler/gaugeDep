args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]

library(qs)
library(dplyr)
library(tidyr)

extract_median_params_ang <- function(dep_type, dep_level, data_num) {
  params <- qread(sprintf("samplers/nimble/sb_mcmc_fits/%s/%s_%s.qs",
                          dep_type, dep_level, data_num)) |> 
    as_tibble() |>
    select(matches("probs|alphastar|betastar")) 
  
  return(params |>
           colMeans() |> 
           t() |>
           as_tibble() |>
           mutate(dataset = data_num))
}


create_tib_med_params_ang <- function(dep_type, dep_level) {
  data_num <- 1:100
  tib_med_params <- sapply(data_num, 
                           function(x) extract_median_params_ang(dep_type = dep_type, 
                                                                 dep_level = dep_level,
                                                                 data_num = x), 
                           simplify = FALSE)
  return(tib_med_params |> bind_rows())
}

med_params <- create_tib_med_params_ang(dep_type = dep_type, 
                                        dep_level = dep_level)

filepath <- sprintf("fits_and_weights/med_params_joint/%s_%s_ang_sb.qs", 
                    dep_type, dep_level)
qsave(med_params, filepath)
print("Posterior means of stick breaking angular parameters have been successfully saved")

