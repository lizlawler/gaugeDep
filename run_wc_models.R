library(geometricMVE)
library(tidyverse)
library(progressr)

load_data <- function(dep_type, dep_level, dataset_num) {
  temp_data <- RcppSimdJson::fload(paste0("data/", dep_type, "/", dep_level, "_", dataset_num, ".json"))
  return(list(r = temp_data$R, w = temp_data$W))
}

grab_exc <- function(r, w) {
  qr <- QR.2d(r = r,w = w, method = "empirical")
  excind <- r > qr$r0w
  rexc <- r[excind]
  wexc <- w[excind]
  
  # Threshold value corresponding to each w:
  r0w <- qr$r0w[excind]
  
  # Check for and remove any NAs - small number can occur because of empirical gauge estimation
  na.ind <- which(is.na(rexc))
  if(length(na.ind) > 0){
    rexc <- rexc[-na.ind]
    wexc <- wexc[-na.ind]
    r0w <- r0w[-na.ind]}
  
  return(list(rexc = rexc, wexc = wexc, r0w=r0w))
}

one_fit <- function(rexc, wexc, r0w, gauge) {
  return(fit.geometricMVE.2d(r = rexc, 
                             w = wexc,
                             r0w = r0w,
                             init.val = c(2, 0.5),
                             gfun = gauge))
}

all_fits <- function(rexc, wexc, r0w) {
  gauge_names <- c("gauge_rvad", "gauge_gaussian", "gauge_invlogistic", "gauge_square")
  temp <- t(sapply(gauge_names, function(x) one_fit(rexc, wexc, r0w, x), USE.NAMES = T)) |>
    as.data.frame() |> 
    rownames_to_column() |> 
    as_tibble()
  best_fit <- temp |> 
    mutate(nllh = unlist(nllh)) |> 
    slice_min(nllh) |> 
    select(rowname, mle)
  return(best_fit)
}

data_and_fit <- function(dep_type, dep_level, dataset_num) {
  data_temp <- load_data(dep_type, dep_level, dataset_num)
  exc_temp <- grab_exc(data_temp$r, data_temp$w)
  return(all_fits(exc_temp$rexc, exc_temp$wexc, exc_temp$r0w))
}

## uncomment the below lines of code to run Wadsworth and Campbell (2023) models --------
# # types <- c("gauss", "logistic")
# types <- "husler_reiss"
# levels <- c("low", "mid", "high")
# datasets <- 1:100
# all_combos <- expand_grid(types, levels, datasets)
# 
# # Enable global progression handlers
# handlers(global = TRUE)
# 
# wc_results <- with_progress({
#   # Create a progress handler
#   p <- progressor(steps = nrow(all_combos))
# 
#   # Apply the function using pmap and update the progress bar
#   pmap(list(all_combos$types, all_combos$levels, all_combos$datasets),
#        ~{
#          p()
#          data_and_fit(..1, ..2, ..3)
#        })
# })

# all_combos <- all_combos |>
#   mutate(gauge_name = map_chr(wc_results, "rowname"),
#          mle = map(wc_results, "mle"))
# saveRDS(all_combos, file = "wadsworth_campbell_hr_fits.RDS")

wc_fits <- readRDS("fits_and_weights/old/wadsworth_campbell_fits.RDS")
wc_hr_fits <- readRDS("fits_and_weights/old/wadsworth_campbell_hr_fits.RDS")
wc_fits_all <- rbind(wc_fits, wc_hr_fits)
split_wc_fits <- split(wc_fits_all, f = list(wc_fits_all$types, wc_fits_all$levels)) |>
  lapply(function(y) y |> mutate(gauge_name = case_when(gauge_name == 'gauge_square' ~ 'rectangular',
                                                        gauge_name == 'gauge_gaussian' ~ 'gauss',
                                                        gauge_name == 'gauge_rvad' ~ 'logistic',
                                                        gauge_name == 'gauge_invlogistic' ~ 'inv_log')))
