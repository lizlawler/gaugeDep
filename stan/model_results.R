library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)

gauge <- "logistic"
dataset <- "high"

csvfiles <- paste0("stan/csv_fits/", gauge, "/",
                   list.files(path = paste0("stan/csv_fits/", gauge, "/"), pattern = paste0(gauge, "_", dataset)))

csvfiles <- csvfiles[!grepl("gauss", csvfiles)]
fit <- as_cmdstan_fit(csvfiles)
fit$diagnostic_summary()

fitmcmc <- as_mcmc.list(fit)

MCMCtrace(fitmcmc,
          params = c('dep1', 'dep2'), 
          gvals = c(0.5, 0.5),
          ind = TRUE, 
          open_pdf = TRUE, 
          filename = paste0("traceplots/", gauge, "/", gauge, "_", dataset, "_alpha.pdf"))
