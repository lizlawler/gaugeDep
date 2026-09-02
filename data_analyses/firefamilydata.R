# =============================================================================
# Preprocessing for the three real fire-weather station datasets (Redstone,
# Thomes Creek, and Friend Mountain RAWS). Reads the raw ERC/FWI CSV exports
# from FireFamily Plus, restricts to the summer fire season (June-October) and
# to nonzero ERC/FWI, and writes each station out as a Stan data list in JSON.
#
# NOTE: the CSV paths below point at a local Desktop folder and must be updated
# to re-run this on another machine. The processed JSON files are already in
# data/raw/, so this only needs re-running if the raw exports change.
#
# Inputs:    local FireFamily Plus CSV exports (see paths below)
# Outputs:   data/raw/erc_fwi_{redstone,thomescreek,friendmtn}.json
#
# The transformation to exponential margins happens later, in
# posterior_marg_transform.R, which writes the *_expo.qs files.
# =============================================================================

library(readr)
library(tidyverse)
library(lubridate)
library(extRemes)

# Shared preprocessing: parse the date, split it into year/month/day, and keep
# only summer-season days with both indices nonzero.
summer_subset <- function(raw) {
  raw |>
    as_tibble() |>
    mutate(date = mdy(DATE)) |>
    separate_wider_delim(cols = date, delim = "-", names = c("yr", "mon", "day")) |>
    mutate(yr  = as.integer(yr),
           mon = as.integer(mon),
           day = as.integer(day)) |>
    filter(mon >= 6, mon <= 10, ERC != 0, FWI != 0) |>
    mutate(date = mdy(DATE))
}

# Package a station's ERC/FWI series as the Stan data list and write to JSON.
write_station_json <- function(station_summer, filename) {
  stan_data <- list(
    N   = as.integer(nrow(station_summer)),
    D   = 2,
    erc = station_summer$ERC,
    fwi = station_summer$FWI
  )
  cmdstanr::write_stan_json(stan_data, sprintf("data/raw/%s", filename))
  return(stan_data)
}

# Quick marginal check: fit a GPD above the 95th percentile and report the
# parameter estimates (used to sanity-check tail behaviour, not saved).
gpd_check <- function(x) {
  thresh95 <- quantile(x, 0.95)
  fit <- fevd(x = x, threshold = thresh95, type = "GP", time.units = "days")
  t(as.matrix(fit$results$par))
}

## Redstone RAWS -- near the Cameron Peak fire --------------------------------
redstone <- read_csv("~/Desktop/firefamily/Redstone_ERC_FWI_take4.csv",
                     na = "NA", comment = "--")
redstone_summer <- summer_subset(redstone)
stan_data_red   <- write_station_json(redstone_summer, "erc_fwi_redstone.json")
gpd_check(stan_data_red$erc)

## Thomes Creek RAWS -- inside the August Complex fire perimeter --------------
thomescreek <- read_csv("~/Desktop/firefamily/040816_ThomesCreek_ERC_FWI.csv",
                        na = "NA", comment = "--")
thomescreek_summer <- summer_subset(thomescreek)
stan_data_tc       <- write_station_json(thomescreek_summer, "erc_fwi_thomescreek.json")
gpd_check(thomescreek_summer$ERC)

## Friend Mountain RAWS -------------------------------------------------------
friendmtn <- read_csv("~/Desktop/firefamily/040512_FriendMtn_ERC_FWI.csv",
                      na = "NA", comment = "--")
friendmtn_summer <- summer_subset(friendmtn)
stan_data_friend <- write_station_json(friendmtn_summer, "erc_fwi_friendmtn.json")
gpd_check(stan_data_friend$fwi)
