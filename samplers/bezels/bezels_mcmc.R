args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]

library(dplyr)
library(tidyr)
library(qs)
library(BezELS)

for(data_num in 1:200) {
  datafile <- sprintf("data/%s/%s_%s.json",
                      dep_type, dep_level, data_num)
  data <- RcppSimdJson::fload(datafile)
  idx <- data$idx
  w <- data$W[idx]
  r <- data$R[idx]
  r0w <- data$r0_w[idx]
  
  results <- fit_mcmc_bezier(N = data$n0,
                             r = r,
                             w = w,
                             r_0 = r0w,
                             iters = 11000, burn = 100, 
                             traceplot = FALSE,
                             print.every = 100)
  qsave(x = results, file = sprintf("samplers/bezels/radial_bezels_fits/%s/%s_%s.qs",
                                    dep_type, dep_level, data_num))
  print(paste0("Successfully saved Bezels fit for dataset number: ", data_num))
}
