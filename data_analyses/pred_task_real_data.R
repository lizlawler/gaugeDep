# =============================================================================
# Creates the prediction task visualisation for the real fire weather data:
# plots the data cloud, threshold, and three prediction boxes (b1, b2, b3)
# representing different joint extremal scenarios.
#
# Inputs:    data/raw/{data_type}_expo.qs
#            data_analyses/egpd_functions.R
# Outputs:   figures/pred_tasks.png (and variants)
# =============================================================================

# library(evd)
# library(extRemes)
library(mvtnorm)
library(tidyverse)
library(patchwork)
library(grafify)
library(qs)
library(cmdstanr)
library(tidyverse)
library(posterior)
source("data_analyses/egpd_functions.R")

# functions to convert back to original data scale 
# CDF of an exponential truncated above at xmax.
ptexp <- function(data, xmax, rate) {
  if(data > xmax) {
    return(1)
  } else {
    numer <- pexp(data, rate)
    denom <- pexp(xmax, rate)
    return(numer/denom)
  }
}

# CDF of the fitted EGPD/exponential mixture used for the FWI margin.
mix_cdf_trunc_expo <- function(data, trunc_pt, params) {
  egpd_part <- g1_cdf(data, sigma = as.numeric(params["sigma"]), xi = as.numeric(params["xi"]), kappa = as.numeric(params["kappa"]))
  expo_part <- rep(NA, length(data))
  for(i in seq_along(data)) {
    expo_part[i] <- ptexp(data[i], trunc_pt, as.numeric(params["rate"]))
  }
  return(as.numeric(params["pi_prob"]) * expo_part + (1 - as.numeric(params["pi_prob"])) * egpd_part)
}

# Root-finding helper: solve mix_cdf_trunc_expo(x) = q for x.
q_mix_cdf <- function(x, q, trunc_pt, params) {
  return(mix_cdf_trunc_expo(x, trunc_pt, params) - q)
}

# Map an exponential-margin value back to the original ERC scale.
erc_og_scale <- function(exp_val, params) {
  return(g1_icdf(pexp(exp_val), 
                 sigma = as.numeric(params["sigma"]), 
                 xi = as.numeric(params["xi"]), 
                 kappa = as.numeric(params["kappa"])))
}

# Map an exponential-margin value back to the original FWI scale, inverting
# the truncated-mixture CDF numerically.
fwi_og_scale <- function(exp_val, params, trunc_pt, f_upper) {
  return(uniroot(q_mix_cdf, 
                 q = pexp(exp_val), 
                 trunc_pt = trunc_pt, 
                 params = params, interval = c(1, f_upper), 
                 extendInt = "up")$root)
}

# friend mountain plots
data_type <- "friendmtn"
friend_data <- qread(sprintf("data/raw/%s_expo.qs", data_type))$cloud_tib
friend_exp_plot <- friend_data |> as_tibble() |> ggplot(aes(x = x, y = y)) + 
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_x_continuous(limits = c(0, 15), expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(limits = c(0, 15), expand = expansion(mult = c(0,0.01))) +
  annotate("rect", xmin = 0.75, xmax = 1.5, ymin = 7,  ymax = 14, fill = get_graf_colours("ok_bluegreen"), color = "black", alpha = 0.75) +
  annotate("rect", xmin = 5.5, xmax = 10.5, ymin = 7,  ymax = 12, fill = get_graf_colours("ok_orange"), color = "black", alpha = 0.75) +
  annotate("rect", xmin = 5, xmax = 11, ymin = 0.25,  ymax = 1, fill = get_graf_colours("ok_redpurple"), color = "black", alpha = 0.75) +
  xlab(expression("ERC")) + ylab(expression("FWI"))
friend_exp_plot

friend_boxes <- list(b1 = list(dim1 = c(0.75, 1.5), dim2 = c(7, 14)),
                     b2 = list(dim1 = c(5.5, 10.5), dim2 = c(7, 12)),
                     b3 = list(dim1 = c(5, 11), dim2 = c(0.25, 1)))

## plot pred task on original scale
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = data_type,
                       full.names = TRUE)
csvfiles <- csvfiles[grepl("trunc45", csvfiles)]
fit <- read_cmdstan_csv(csvfiles, variables = c("xi", "kappa", "sigma", "rate", "pi_prob"))$post_warmup_draws

params <- fit |> as_draws_df() |>
  rename(draw = ".draw") |>
  select(!contains(c("log", "lp", "chain", "iter"))) |>
  pivot_longer(cols = -"draw") |>
  separate_wider_delim(cols = "name", delim = "[", names = c("param", "index"), too_few = "align_start") |>
  mutate(index = case_when(grepl("1", index) ~ "erc",
                           grepl("2", index) ~ "fwi",
                           is.na(index) ~ "fwi"),
         value = case_when(param == "pi_prob" ~ (1 / (1 + exp(-value))),
                           .default = value))

friend_erc_params <- params |> filter(index == "erc") |> pivot_wider(names_from = "param", values_from = "value") |> select(-c(index, draw)) |> colMeans()
friend_fwi_params <- params |> filter(index == "fwi") |> pivot_wider(names_from = "param", values_from = "value") |> select(-c(index, draw)) |> colMeans()

og_data <- RcppSimdJson::fload(sprintf("data/erc_fwi_%s.json", data_type))
friend_og_tib <- tibble(erc = og_data$erc, fwi = og_data$fwi)

trunc_friend <- 45
upper_friend <- 65
friend_og_plot <- friend_og_tib |> ggplot(aes(x = erc, y = fwi)) + 
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_x_continuous(expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  annotate("rect", 
           xmin = erc_og_scale(0.75, friend_erc_params), xmax = erc_og_scale(1.5, friend_erc_params), 
           ymin = fwi_og_scale(7, friend_fwi_params, trunc_friend, upper_friend),  ymax = fwi_og_scale(14, friend_fwi_params, trunc_friend, upper_friend), 
           fill = get_graf_colours("ok_bluegreen"), color = "black", alpha = 0.75) +
  annotate("rect", 
           xmin = erc_og_scale(5.5, friend_erc_params), xmax = erc_og_scale(10.5, friend_erc_params), 
           ymin = fwi_og_scale(7, friend_fwi_params, trunc_friend, upper_friend),  ymax = fwi_og_scale(12, friend_fwi_params, trunc_friend, upper_friend), 
           fill = get_graf_colours("ok_orange"), color = "black", alpha = 0.75) +
  annotate("rect", 
           xmin = erc_og_scale(5, friend_erc_params), xmax = erc_og_scale(11, friend_erc_params), 
           ymin = fwi_og_scale(0.25, friend_fwi_params, trunc_friend, upper_friend),  ymax = fwi_og_scale(1, friend_fwi_params, trunc_friend, upper_friend), 
           fill = get_graf_colours("ok_redpurple"), color = "black", alpha = 0.75) +
  xlab(expression("ERC")) + ylab(expression("FWI"))
friend_og_plot
friend_plots <- (friend_og_plot + plot_spacer() + friend_exp_plot + plot_layout(widths = c(0.95, 0.025, 0.95))) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
friend_plots

## redstone plots
data_type <- "redstone"
redstone_data <- qread(sprintf("data/raw/%s_expo.qs", data_type))$cloud_tib

redstone_exp_plot <- redstone_data |> as_tibble() |> ggplot(aes(x = x, y = y)) + 
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_x_continuous(limits = c(0,15), expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(limits = c(0, 15), expand = expansion(mult = c(0,0.01))) +
  annotate("rect", xmin = 0.75, xmax = 1.5, ymin = 9,  ymax = 14, fill = get_graf_colours("ok_bluegreen"), color = "black", alpha = 0.75) +
  annotate("rect", xmin = 6, xmax = 11, ymin = 7,  ymax = 12, fill = get_graf_colours("ok_orange"), color = "black", alpha = 0.75) +
  annotate("rect", xmin = 6, xmax = 12, ymin = 0.25,  ymax = 1, fill = get_graf_colours("ok_redpurple"), color = "black", alpha = 0.75) +
  xlab(expression("ERC")) + ylab(expression("FWI"))
redstone_exp_plot

redstone_boxes <- list(b1 = list(dim1 = c(0.75, 1.5), dim2 = c(9, 14)),
                       b2 = list(dim1 = c(6, 11), dim2 = c(7, 12)),
                       b3 = list(dim1 = c(6, 12), dim2 = c(0.25, 1)))

## plot pred task on original scale
csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
                       pattern = data_type,
                       full.names = TRUE)
csvfiles <- csvfiles[grepl("v2", csvfiles)]
fit <- read_cmdstan_csv(csvfiles, variables = c("xi", "kappa", "sigma", "rate", "pi_prob"))$post_warmup_draws

params <- fit |> as_draws_df() |>
  rename(draw = ".draw") |>
  select(!contains(c("log", "lp", "chain", "iter"))) |>
  pivot_longer(cols = -"draw") |>
  separate_wider_delim(cols = "name", delim = "[", names = c("param", "index"), too_few = "align_start") |>
  mutate(index = case_when(grepl("1", index) ~ "erc",
                           grepl("2", index) ~ "fwi",
                           is.na(index) ~ "fwi"),
         value = case_when(param == "pi_prob" ~ (1 / (1 + exp(-value))),
                           .default = value))

redstone_erc_params <- params |> filter(index == "erc") |> pivot_wider(names_from = "param", values_from = "value") |> select(-c(index, draw)) |> colMeans()
redstone_fwi_params <- params |> filter(index == "fwi") |> pivot_wider(names_from = "param", values_from = "value") |> select(-c(index, draw)) |> colMeans()

og_data <- RcppSimdJson::fload(sprintf("data/erc_fwi_%s.json", data_type))
redstone_og_tib <- tibble(erc = og_data$erc, fwi = og_data$fwi)

trunc_red <- 100
upper_red <- 120
redstone_og_plot <- redstone_og_tib |> ggplot(aes(x = erc, y = fwi)) + 
  geom_point() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent', color = 'transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_x_continuous(expand = expansion(mult = c(0,0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.01))) +
  annotate("rect", 
           xmin = erc_og_scale(0.75, redstone_erc_params), xmax = erc_og_scale(1.5, redstone_erc_params), 
           ymin = fwi_og_scale(9, redstone_fwi_params, trunc_red, upper_red),  ymax = fwi_og_scale(14, redstone_fwi_params, trunc_red, upper_red), 
           fill = get_graf_colours("ok_bluegreen"), color = "black", alpha = 0.75) +
  annotate("rect", 
           xmin = erc_og_scale(6, redstone_erc_params), xmax = erc_og_scale(11, redstone_erc_params), 
           ymin = fwi_og_scale(7, redstone_fwi_params, trunc_red, upper_red),  ymax = fwi_og_scale(12, redstone_fwi_params, trunc_red, upper_red), 
           fill = get_graf_colours("ok_orange"), color = "black", alpha = 0.75) +
  annotate("rect", 
           xmin = erc_og_scale(6, redstone_erc_params), xmax = erc_og_scale(12, redstone_erc_params), 
           ymin = fwi_og_scale(0.25, redstone_fwi_params, trunc_red, upper_red),  ymax = fwi_og_scale(1, redstone_fwi_params, trunc_red, upper_red), 
           fill = get_graf_colours("ok_redpurple"), color = "black", alpha = 0.75) +
  xlab(expression("ERC")) + ylab(expression("FWI")) 
redstone_og_plot
redstone_plots <- (redstone_og_plot + plot_spacer() + redstone_exp_plot + plot_layout(widths = c(0.95, 0.025, 0.95))) & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
redstone_plots


friend_plots / plot_spacer() / redstone_plots + plot_layout(heights = c(0.95,0.025, 0.95))  & 
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color='transparent'))
ggsave("figures/both_pred_tasks.png",
       dpi = 320,
       width = 10,
       height = 10,
       bg = 'transparent')
knitr::plot_crop("figures/both_pred_tasks.png")

all_boxes <- list(friendmtn = friend_boxes, redstone = redstone_boxes)
qsave(x = all_boxes, "data_analyses/pred_boxes.qs")
