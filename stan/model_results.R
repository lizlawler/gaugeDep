library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)

gauge <- "gauss"
dataset <- "high"

csvfiles <- paste0("stan/csv_fits/calibrate/", gauge, "/",
                   list.files(path = paste0("stan/csv_fits/calibrate/", gauge, "/"), 
                              pattern = paste0(dataset, "_10_cens_marg_\\d{1}.csv")))

fit <- as_cmdstan_fit(csvfiles)
fit$diagnostic_summary()

fitmcmc <- as_mcmc.list(fit)

MCMCtrace(fitmcmc,
          params = c('alpha', 'dep'), 
          gvals = c(2, .9),
          ind = TRUE, 
          open_pdf = TRUE)
          # filename = paste0("traceplots/", gauge, "/", gauge, "_", dataset, "_alpha.pdf"))
