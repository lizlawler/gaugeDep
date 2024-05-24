args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
gauge <- args[2]
likelihood <- args[3]
threshold <- args[4]
dep_level <- args[5]

library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)

extract_median_params <- function(sim_phase = "stacking", gauge, dep_type, likelihood, threshold, dep_level, dataset_num) {
  start_file_path <- paste0("stan/csv_fits/", sim_phase, "/", dep_type, "/", gauge, "/")
  csvfiles <- paste0(start_file_path,
                     list.files(path = start_file_path, 
                                pattern = paste0(dep_level, "_", dataset_num, "_", likelihood, "_", threshold, "_\\d{1}.csv")))
  if(gauge != "dirichlet") {
    return(fit |> as_draws_df() |> 
             select(alpha, dep) |> 
             apply(MARGIN = 2, FUN = median) |> t() |>
             as_tibble() |>
             mutate(dataset = dataset_num))
  } else {
    return(fit |> as_draws_df() |> 
             select(alpha, theta1, theta2) |> 
             apply(MARGIN = 2, FUN = median) |> t() |>
             as_tibble() |>
             mutate(dataset = dataset_num))
  }
}

create_tib_med_params <- function(sim_phase = "stacking", gauge, dep_type, likelihood, threshold, dep_level) {
  dataset_num <- 1:100
  tib_med_params <- sapply(dataset_num, 
                           function(x) extract_median_params(gauge = gauge, 
                                                             dep_type = dep_type,
                                                             likelihood = likelihood, 
                                                             threshold = threshold, 
                                                             dep_level = dep_level, 
                                                             dataset_num = x), 
                           simplify = FALSE)
  return(tib_med_params |> bind_rows())
}

med_params <- create_tib_med_params("stacking", gauge = gauge, dep_type = dep_type,
                                    likelihood = likelihood, threshold = threshold, dep_level = dep_level)

filepath <- paste0("extracted_params/", gauge, "_", dep_type, "_", dep_level, "_", likelihood, "_", threshold, "_params.RDS")
saveRDS(med_params, filepath)

print("Posterior medians of parameters have been successfully saved")
