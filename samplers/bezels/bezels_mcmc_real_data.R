library(dplyr)
library(tidyr)
library(qs)
library(BezELS)

data_type <- "redstone"

data <- qread(sprintf("data/%s_expo.qs", data_type))
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
qsave(x = results, file = sprintf("samplers/bezels/radial_bezels_fits/real_data/%s.qs", data_type))
print(paste0("Successfully saved Bezels fit for on fire data"))
