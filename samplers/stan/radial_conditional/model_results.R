library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)

create_mcmc_object <- function(gauge, dep_level, lhood_type, thresh_type, data_num) {
  csvfiles <- paste0("stan/csv_fits/calibrate/", gauge, "/",
                     list.files(path = paste0("stan/csv_fits/calibrate/", gauge, "/"), 
                                pattern = paste0(dep_level, "_", data_num, "_", lhood_type, "_", thresh_type, "_\\d{1}.csv")))
  fit <- as_cmdstan_fit(csvfiles)
  return(as_mcmc.list(fit))
}

save_mcmc_trace <- function(gauge, dep_level, lhood_type, thresh_type, data_num, truth) {
  temp <- create_mcmc_object(gauge, dep_level, lhood_type, thresh_type, data_num)
  if (gauge != "dirichlet") {
    MCMCtrace(temp,
              params = c('alpha', 'dep'), 
              gvals = c(2, truth),
              ind = TRUE, 
              open_pdf = FALSE,
              filename = paste0("traceplots/", gauge, "/",dep_level, "_", lhood_type, "_", thresh_type, "_", data_num, ".pdf"))
    print(paste0(data_num, " of 100 complete"))
  } else {
    MCMCtrace(temp,
              params = c('alpha', 'theta1', 'theta2'), 
              gvals = c(2, truth, 2),
              ind = TRUE, 
              open_pdf = FALSE,
              filename = paste0("traceplots/", gauge, "/", dep_level, "_", lhood_type, "_", thresh_type, "_", data_num, ".pdf"))
    print(paste0(data_num, " of 100 complete"))
  }
}

true_vals <- c(0.1, 0.5, 0.9)
dep_names <- c("low", "mid", "high")
for(i in seq_along(true_vals)) {
  sapply(1:100, function(y) save_mcmc_trace("gauss", dep_names[i], "cens", "marg_take2", y, true_vals[i]))
}
