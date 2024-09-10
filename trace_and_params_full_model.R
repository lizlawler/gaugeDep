args <- commandArgs(trailingOnly=TRUE)
gauge <- args[1]
dep_level <- args[2]

library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)

filepath <- paste0("stan/radial_angular/csv_fits/", gauge, "/")
sapply(1:100, function(i) {
  files <- list.files(path = filepath, 
                      pattern = paste0(dep_level, "_", i, "_\\d{1}.csv"), full.names = TRUE)
  fit <- as_cmdstan_fit(files)
  MCMCtrace(as_mcmc.list(fit),
            ind = TRUE, 
            open_pdf = FALSE,
            filename = paste0("trace_and_params_mix_betas/", gauge, "/traceplots/", dep_level, "_", i, ".pdf"))
  med_params <- fit |> as_draws_df() |> 
    select(any_of(contains(c("weights", "alpha","beta", "dep", ".draw", ".chain")))) |>
    group_by(.chain) |> 
    summarize(across(everything(), median)) |> 
    select(-.draw)
  saveRDS(med_params, file = paste0("trace_and_params_mix_betas/", gauge, "/params/", dep_level, "_", i, ".rds"))
})

print(paste0("Traceplots and posterior parameter estimates have been saved to disk for ", gauge, ", ", dep_level))