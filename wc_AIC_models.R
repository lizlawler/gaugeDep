library(tidyverse)
wadsworth_campbell_fits <- readRDS("~/Desktop/research/gaugeDependence/wadsworth_campbell_fits.RDS")

fits <- split(wadsworth_campbell_fits, ~ types + levels)
table(fits[["gauss.high"]]$gauge_name)
table(fits[["gauss.mid"]]$gauge_name)
table(fits[["gauss.low"]]$gauge_name)

table(fits[["logistic.high"]]$gauge_name)
table(fits[["logistic.mid"]]$gauge_name)
table(fits[["logistic.low"]]$gauge_name)


wc_hr_fits <- readRDS("~/Desktop/research/gaugeDependence/wadsworth_campbell_hr_fits.RDS")
hr_fits <- split(wc_hr_fits, ~ types + levels)
table(hr_fits[["husler_reiss.high"]]$gauge_name)
table(hr_fits[["husler_reiss.mid"]]$gauge_name)
table(hr_fits[["husler_reiss.low"]]$gauge_name)
