args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]

library(qs)
library(dplyr)
library(tidyr)

extract_post_params_ang_mix <- function(dep_type, dep_level, data_num) {
  params <- qread(sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
                          dep_type, dep_level, data_num)) |>
    as_tibble() |>
    select(matches("probs|alphastar|betastar"))
  
  return(params |>
           colMeans() |>
           t() |>
           as_tibble() |>
           mutate(dataset = data_num))
}


create_tib_post_params_ang_mix <- function(dep_type, dep_level) {
  data_num <- 1:200
  tib_post_params <- sapply(data_num,
                            function(x) extract_post_params_ang_mix(dep_type = dep_type,
                                                                    dep_level = dep_level,
                                                                    data_num = x),
                            simplify = FALSE) |> 
    bind_rows()
  filepath <- sprintf("fits_and_weights/post_params_joint/%s_%s_ang_mix.qs",
                      dep_type, dep_level)
  qsave(tib_post_params, filepath)
  print("Posterior means of stick breaking angular parameters have been successfully saved")
}

create_tib_post_params_ang_mix(dep_type = dep_type,
                               dep_level = dep_level)



