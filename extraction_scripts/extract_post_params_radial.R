args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
gauge <- args[3]
likelihood <- args[4]

library(qs)
library(dplyr)
library(tidyr)

extract_post_params_radial <- function(dep_type, dep_level, gauge, data_num, likelihood) {
  params <- qread(sprintf("samplers/rcpp/radial_mcmc_fits/%s/%s_%s_%s_%s.qs",
                          dep_type, gauge, likelihood, dep_level, data_num))$samples |> as_tibble()
  return(params |> 
           select(-sigma_m) |>
           apply(MARGIN = 2, FUN = median) |> 
           t() |>
           as_tibble() |>
           mutate(dataset = data_num))
}


create_tib_post_params_radial <- function(dep_type, dep_level, gauge, likelihood) {
  data_num <- 1:200
  tib_post_params <- sapply(data_num, 
                            function(x) extract_post_params_radial(dep_type = dep_type, 
                                                                   dep_level = dep_level, 
                                                                   gauge = gauge, 
                                                                   data_num = x, 
                                                                   likelihood = likelihood), 
                            simplify = FALSE) |> 
    bind_rows()
  filepath <- sprintf("fits_and_weights/post_params_joint/%s_%s_%s_%s_radial.qs", 
                      gauge, dep_type, dep_level, likelihood)
  qsave(tib_post_params, filepath)
  print("Posterior medians of radial parameters have been successfully saved")
}

create_tib_post_params_radial(dep_type = dep_type, 
                              dep_level = dep_level, 
                              gauge = gauge, 
                              likelihood = likelihood)
