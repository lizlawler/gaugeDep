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

print(paste0("dep_type = ", dep_type))
print(paste0("gauge = ", gauge))
print(paste0("likelihood = ", likelihood))
print(paste0("threshold = ", threshold))
print(paste0("dep_level = ", dep_level))

extract_params <- function(sim_phase = "stacking", gauge, dep_type, likelihood, threshold, dep_level, dataset_num) {
  start_file_path <- paste0("stan/csv_fits/", sim_phase, "/", dep_type, "/", gauge, "/")
  csvfiles <- paste0(start_file_path,
                     list.files(path = start_file_path, 
                                pattern = paste0(dep_level, "_", dataset_num, "_", likelihood, "_", threshold, "_\\d{1}.csv")))
  if(gauge != "dirichlet") {
    fit <- read_cmdstan_csv(csvfiles, variables = c("alpha", "dep"))$post_warmup_draws
    return(fit |> as_draws_df() |> 
             rename(draw = ".draw") |>
             select(alpha, dep, draw) |>
             as_tibble() |>
             mutate(dataset = dataset_num))
  } else {
    fit <- read_cmdstan_csv(csvfiles, variables = c("alpha", "theta1", "theta2"))$post_warmup_draws
    return(fit |> as_draws_df() |> 
             rename(draw = ".draw") |>
             select(alpha, theta1, theta2, draw) |>
             as_tibble() |>
             mutate(dataset = dataset_num))
  }
}

create_tib_params <- function(sim_phase = "stacking", gauge, dep_type, likelihood, threshold, dep_level) {
  data_num <- 1:100
  tib_params <- sapply(data_num, 
                       function(x) extract_params(gauge = gauge, 
                                                  dep_type = dep_type,
                                                  likelihood = likelihood, 
                                                  threshold = threshold, 
                                                  dep_level = dep_level, 
                                                  dataset_num = x), 
                       simplify = FALSE)
  return(tib_params)
}

all_iter_params <- create_tib_params("stacking", gauge = gauge, dep_type = dep_type,
                                     likelihood = likelihood, threshold = threshold, dep_level = dep_level)

filepath <- paste0("extracted_params/", gauge, "_", dep_type, "_", dep_level, "_", likelihood, "_", threshold, "_all_iter_params.RDS")
saveRDS(all_iter_params, filepath)

print("Posterior samples of parameters have been successfully saved")
