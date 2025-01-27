library(readr)
library(tidyverse)
library(lubridate)
# library(cffdrs)
# redstone_weather <- read_csv("~/Desktop/firefamily/firefamily/050508_Redstone_weather.csv", col_names = FALSE)
# colnames(redstone_weather) <- c("id", "date", "time", "temp", "rh", "prec", "ws")
# redstone_weather <- redstone_weather |> as_tibble() |> 
#   mutate(temp = (temp - 32) * (5/9),
#          ws = 1.609344 * ws,
#          date = mdy(date)) |>
#   separate_wider_delim(cols = date, delim = "-", names = c("yr", "mon", "day")) |>
#   mutate(yr = as.integer(yr),
#          mon = as.integer(mon),
#          day = as.integer(day),
#          lat = 40.57)
# # redstone_weather$rh[4022] <- 99
# # fwi_redstone <- fwi(redstone_weather)

redstone <- read_csv("~/Desktop/firefamily/data/050508_Redstone_daily_ERC_FWI.csv",
                       na = "NA", comment = "--")
redstone_nonfire <- redstone |> as_tibble() |> 
  mutate(date = mdy(DATE)) |>
  separate_wider_delim(cols = date, delim = "-", names = c("yr", "mon", "day")) |>
  mutate(yr = as.integer(yr),
         mon = as.integer(mon),
         day = as.integer(day)) |>
  filter(mon < 6 | mon > 10)

redstone_summer <- redstone |> as_tibble() |> 
  mutate(date = mdy(DATE)) |>
  separate_wider_delim(cols = date, delim = "-", names = c("yr", "mon", "day")) |>
  mutate(yr = as.integer(yr),
         mon = as.integer(mon),
         day = as.integer(day)) |>
  filter(mon >= 6, mon <= 10) |>
  mutate(date = mdy(DATE))
## make stan data list
# erc_fwi <- cbind(redstone_summer$ERC, redstone_summer$FWI) |> t()
idx <- unique(c(which(erc_fwi[1,] == 0), which(erc_fwi[2,] == 0)))
plot(redstone_summer$ERC[-idx], redstone_summer$FWI[-idx])
N <- as.integer(nrow(redstone_summer) - length(idx))
erc <- redstone_summer$ERC[-idx]
fwi <- redstone_summer$FWI[-idx]
plot(erc, fwi)

stan_data_og <- list(
  N = N,
  D = 2,
  erc = erc,
  fwi = fwi
)
cmdstanr::write_stan_json(stan_data_og, "data/erc_fwi_og.json")

stan_data_std <- list(
  N = N,
  D = 2,
  erc = erc / sd(erc),
  fwi = fwi / sd(fwi)
)
cmdstanr::write_stan_json(stan_data_std, "data/erc_fwi_std.json")
# ## rank transform to expo margins
exp_ERC_rank <- qexp(rank(erc_data)/(length(erc_data) + 1))
exp_FWI_rank <- qexp(rank(fwi_data)/(length(fwi_data) + 1))
# 
# exp_ERC_nofire <- qexp(rank(redstone_nonfire$ERC)/(nrow(redstone_nonfire) + 1))
# exp_FWI_nofire <- qexp(rank(redstone_nonfire$FWI)/(nrow(redstone_nonfire) + 1))
# 
# plot(exp_ERC_nofire/log(nrow(redstone_nonfire)), exp_FWI_nofire/log(nrow(redstone_nonfire)), pch = 20)
# idx <- as.integer(redstone_summer$date)
# lm(redstone_summer$ERC ~ idx)

