# =============================================================================
# Fits the EGPD marginal model via Stan and plots the posterior predictive
# marginal density alongside the observed ERC/FWI histograms. Also produces
# PP/QQ plots and the quantile assumption diagnostic figure.
#
# Inputs:    samplers/stan/marg_transform/csv_fits/  (Stan output)
#            data/raw/erc_fwi_{station}.json
# Outputs:   figures/marg_transform.pdf, figures/quant_assumpt_plot.pdf
# =============================================================================

library(cmdstanr)
library(tidyverse)
library(MCMCvis)
library(posterior)
library(patchwork)
library(qqplotr)
source("data_analyses/egpd_functions.R")

## Mixture of EGPD with truncated exponential
# Density of an exponential truncated above at xmax.
dtexp <- function(x, xmax, rate) {
  if (x > xmax) {
    return(0)
  } else {
    numer <- dexp(x, rate)
    denom <- pexp(xmax, rate)
    return(numer / denom)
  }
}

# CDF of an exponential truncated above at xmax.
ptexp <- function(data, xmax, rate) {
  if (data > xmax) {
    return(1)
  } else {
    numer <- pexp(data, rate)
    denom <- pexp(xmax, rate)
    return(numer / denom)
  }
}

# Mixture density: EGPD in the upper tail, truncated exponential below the
# threshold. This gives a flexible bulk without imposing a hard threshold.
mix_dens_trunc_expo <- function(data, trunc_pt, params) {
  egpd_part <- g1_pdf(data, sigma = params["sigma"], xi = params["xi"], kappa = params["kappa"])
  expo_part <- rep(NA, length(data))
  for (i in seq_along(data)) {
    expo_part[i] <- dtexp(data[i], trunc_pt, params["rate"])
  }
  return(params["pi_prob"] * expo_part + (1 - params["pi_prob"]) * egpd_part)
}

# CDF of the EGPD / truncated-exponential mixture; used for the probability
# integral transform onto uniform (then exponential) margins.
mix_cdf_trunc_expo <- function(data, trunc_pt, params) {
  egpd_part <- g1_cdf(data, sigma = params["sigma"], xi = params["xi"], kappa = params["kappa"])
  expo_part <- rep(NA, length(data))
  for (i in seq_along(data)) {
    expo_part[i] <- ptexp(data[i], trunc_pt, params["rate"])
  }
  return(params["pi_prob"] * expo_part + (1 - params["pi_prob"]) * egpd_part)
}


# Find the marginal threshold leaving n0 points above it and return the
# per-angle radial threshold r0(w) along with the exceedance indices.
grab_top_n <- function(cloud_tib, n0 = 1, N = 5000) {
  tau <- (N - n0) / N
  # n0 <- ceiling((1 - tau) * N)
  q1 <- quantile(cloud_tib$x, tau)
  q2 <- quantile(cloud_tib$y, tau)
  q <- max(q1, q2)
  idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  eps <- 0.001
  while (length(idx) > n0) {
    q <- q + eps
    idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  }
  cloud_tib <- cloud_tib |>
    mutate(
      r0_w = ifelse(w1 > 0.5, q / w1, q / w2),
      x_lb = ifelse(w1 < 0.5, q, q * y / x),
      y_lb = ifelse(w1 > 0.5, q, q * x / y),
      high = as.factor(case_when(y > q | x > q ~ 1,
        .default = 0
      ))
    )
  return(list(
    q = q,
    idx = idx,
    n0 = length(idx),
    N = N,
    R = cloud_tib$r,
    W = cloud_tib$w1,
    W2 = cloud_tib$w2,
    r0_w = cloud_tib$r0_w,
    cloud_tib = cloud_tib
  ))
}

cols <- c("lightblue", "blue", "red")

# Each station block below follows the same three steps: read the Stan fit,
# transform ERC/FWI to exponential margins via the fitted mixture CDF, and
# save the polar (R, W) representation plus threshold for the MCMC samplers.
### Redstone weather station
## read in stan results
csvfiles <- list.files(
  path = "samplers/stan/marg_transform/csv_fits/",
  pattern = "redstone_transform_v2",
  full.names = TRUE
)
red_fit <- read_cmdstan_csv(csvfiles, variables = c("xi", "kappa", "sigma", "rate", "pi_prob"))$post_warmup_draws

red_params <- red_fit |>
  as_draws_df() |>
  rename(draw = ".draw") |>
  select(!contains(c("log", "lp", "chain", "iter"))) |>
  pivot_longer(cols = -"draw") |>
  separate_wider_delim(cols = "name", delim = "[", names = c("param", "index"), too_few = "align_start") |>
  mutate(
    index = case_when(
      grepl("1", index) ~ "erc",
      grepl("2", index) ~ "fwi",
      is.na(index) ~ "fwi"
    ),
    value = case_when(param == "pi_prob" ~ (1 / (1 + exp(-value))),
      .default = value
    )
  )

erc_params <- red_params |>
  filter(index == "erc") |>
  pivot_wider(names_from = "param", values_from = "value") |>
  select(-c(index, draw)) |>
  colMeans()
fwi_params <- red_params |>
  filter(index == "fwi") |>
  pivot_wider(names_from = "param", values_from = "value") |>
  select(-c(index, draw)) |>
  colMeans()

data <- RcppSimdJson::fload("data/erc_fwi_redstone.json")

hist(data$erc, freq = FALSE, breaks = 35)
curve(g1_pdf(x, sigma = erc_params["sigma"], xi = erc_params["xi"], kappa = erc_params["kappa"]), add = TRUE)
erc_unif <- g1_cdf(data$erc, sigma = erc_params["sigma"], xi = erc_params["xi"], kappa = erc_params["kappa"])
hist(erc_unif)

hist(data$fwi, freq = FALSE, breaks = 45)
curve(mix_dens_trunc_expo(x, 100, fwi_params), add = TRUE)

fwi_unif <- mix_cdf_trunc_expo(data$fwi, 100, fwi_params)
hist(fwi_unif)

fwi_exp <- qexp(fwi_unif)
erc_exp <- qexp(erc_unif)


indices_expo_red <- cbind(x = erc_exp, y = fwi_exp) |>
  as_tibble() |>
  mutate(
    r = x + y,
    w1 = x / r,
    w2 = y / r,
    dataset = "redstone"
  )

fire_data_list <- grab_top_n(indices_expo, n0 = 155, N = nrow(indices_expo))
n <- nrow(fire_data_list$cloud_tib)
fire_data_list$cloud_tib |>
  ggplot(aes(x / log(n), y / log(n), color = high)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(x = y_lb / log(n), y = x_lb / log(n), color = "red"), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  # scale_x_continuous(expand = c(0,0.01)) +
  # scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.025))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  theme(
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.position = "none"
  ) +
  xlab(expression("ERC" / "log(n)")) +
  ylab(expression("FWI" / "log(n)"))

fire_data_list$cloud_tib |>
  ggplot(aes(w1, r / log(n), color = high)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(x = x_lb / (x_lb + y_lb), y = (x_lb + y_lb) / log(n), color = "red"), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  theme(
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.position = "none"
  ) +
  xlab(expression("W")) +
  ylab(expression("R" / "log(n)"))

qs::qsave(fire_data_list, file = "data/raw/redstone_expo.qs")

## plots of transformation
erc_tibble <- tibble(idx_val = data$erc, idx = "erc") |>
  mutate(
    grid_val = seq(min(idx_val), max(idx_val) + 1, length.out = n()),
    dens_val = g1_pdf(grid_val, sigma = erc_params["sigma"], xi = erc_params["xi"], kappa = erc_params["kappa"])
  )

fwi_tibble <- tibble(idx_val = data$fwi, idx = "fwi") |>
  mutate(
    grid_val = seq(min(idx_val), max(idx_val) + 1, length.out = n()),
    dens_val = mix_dens_trunc_expo(grid_val, 100, fwi_params)
  )

red_tibble <- rbind(erc_tibble, fwi_tibble)

erc_plot <- erc_tibble |> ggplot(aes(x = idx_val)) +
  geom_histogram(aes(y = after_stat(density)),
    bins = 25,
    boundary = 0, alpha = 0.5, position = "identity", color = "grey34", fill = "grey"
  ) +
  geom_line(aes(x = grid_val, y = dens_val), color = "blue", linewidth = 1) +
  xlab("ERC") +
  ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "none",
    legend.background = element_rect(fill = "transparent", color = "transparent")
  )

fwi_plot <- fwi_tibble |> ggplot(aes(x = idx_val)) +
  geom_histogram(aes(y = after_stat(density)),
    bins = 30,
    boundary = 0, alpha = 0.5, position = "identity", color = "grey34", fill = "grey"
  ) +
  geom_line(aes(x = grid_val, y = dens_val), color = "blue", linewidth = 1) +
  xlab("FWI") +
  ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "none",
    legend.background = element_rect(fill = "transparent", color = "transparent")
  )

redstone_transforms <- erc_plot + plot_spacer() + fwi_plot + plot_layout(widths = c(0.95, 0.05, 0.95))
# ggsave("figures/fire_hists.pdf",
#        dpi = 320,
#        bg = "transparent",
#        width = 8,
#        height = 4)
# knitr::plot_crop("figures/fire_hists.pdf")

redstone_qq <- indices_expo_red |>
  select(x, y) |>
  rename(ERC = x, FWI = y) |>
  pivot_longer(cols = 1:2, names_to = "idx", values_to = "value") |>
  ggplot(aes(sample = value)) +
  stat_qq_band(distribution = "exp", dparams = list(rate = 1), bandType = "ell", fill = "transparent", color = "darkgrey", linewidth = 0.75) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, col = "red", linetype = "dashed") +
  stat_qq_point(distribution = "exp", dparams = list(rate = 1), alpha = 0.75) +
  scale_x_continuous(limits = c(NA, 13), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(NA, 13), expand = expansion(mult = c(0, 0.01))) +
  facet_grid(~idx, axes = "all", axis.labels = "margins") +
  labs(
    x = "Theoretical exponential quantiles",
    y = "Sample quantiles"
  ) +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.75, "cm"),
    panel.spacing.y = unit(0.75, "cm"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text.x = element_text(size = rel(1.5)),
    strip.text.y = element_blank(),
    axis.text = element_text(size = rel(1.3)),
    axis.title = element_text(size = rel(1.3))
  )

ggsave("figures/qq_plot_marg_transform_redstone.pdf",
  dpi = 320,
  bg = "transparent",
  width = 8,
  height = 4
)
knitr::plot_crop("figures/qq_plot_marg_transform_redstone.pdf")


# ### Thomes Creek weather station
# ## read in stan results
# csvfiles <- list.files(path = "samplers/stan/marg_transform/csv_fits/",
#                        pattern = "thomescreek_g2",
#                        full.names = TRUE)
# tc_fit <- read_cmdstan_csv(csvfiles, variables = c("xi", "kappa1", "kappa2", "sigma", "prob"))$post_warmup_draws
# tc_params <- tc_fit |> as_draws_df() |>
#   rename(draw = ".draw") |>
#   select(!contains(c("log", "lp", "iter"))) |>
#   pivot_longer(cols = -c("draw", ".chain")) |>
#   separate_wider_delim(cols = "name", delim = "[", names = c("param", "index"), too_few = "align_start") |>
#   mutate(index = case_when(grepl("1", index) ~ "erc",
#                            grepl("2", index) ~ "fwi",
#                            is.na(index) ~ "fwi"),
#          value = case_when(param == "pi_prob" ~ (1 / (1 + exp(-value))),
#                            .default = value))
#
# erc_params <- tc_params |> filter(index == "erc") |> select(-index) |>
#   group_by(param, .chain) |> summarize(mean_param = mean(value)) |> ungroup() |> filter(.chain != 1)
# erc_params <- tc_params |> filter(index == "erc", .chain != 1) |> pivot_wider(names_from = "param", values_from = "value") |> select(-c(index, draw, .chain)) |> colMeans()
# fwi_params <- tc_params |> filter(index == "fwi") |> select(-index) |>
#   group_by(param, .chain) |> summarize(mean_param = mean(value)) |> ungroup() |> filter(.chain != 3)
# fwi_params <- tc_params |> filter(index == "fwi", .chain != 3) |> pivot_wider(names_from = "param", values_from = "value") |> select(-c(index, draw, .chain)) |> colMeans()


### Friend Mountain weather station
## read in stan results
csvfiles <- list.files(
  path = "samplers/stan/marg_transform/csv_fits/",
  pattern = "friendmtn_trunc45",
  full.names = TRUE
)
fm_fit <- read_cmdstan_csv(csvfiles, variables = c("xi", "kappa", "sigma", "rate", "pi_prob"))$post_warmup_draws

fm_params <- fm_fit |>
  as_draws_df() |>
  rename(draw = ".draw") |>
  select(!contains(c("log", "lp", "chain", "iter"))) |>
  pivot_longer(cols = -"draw") |>
  separate_wider_delim(cols = "name", delim = "[", names = c("param", "index"), too_few = "align_start") |>
  mutate(
    index = case_when(
      grepl("1", index) ~ "erc",
      grepl("2", index) ~ "fwi",
      is.na(index) ~ "fwi"
    ),
    value = case_when(param == "pi_prob" ~ (1 / (1 + exp(-value))),
      .default = value
    )
  )

erc_params <- fm_params |>
  filter(index == "erc") |>
  pivot_wider(names_from = "param", values_from = "value") |>
  select(-c(index, draw)) |>
  colMeans()
fwi_params <- fm_params |>
  filter(index == "fwi") |>
  pivot_wider(names_from = "param", values_from = "value") |>
  select(-c(index, draw)) |>
  colMeans()

data <- RcppSimdJson::fload("data/erc_fwi_friendmtn.json")

hist(data$erc, freq = FALSE, breaks = 30)
curve(g1_pdf(x, sigma = erc_params["sigma"], xi = erc_params["xi"], kappa = erc_params["kappa"]), add = TRUE)
erc_unif <- g1_cdf(data$erc, sigma = erc_params["sigma"], xi = erc_params["xi"], kappa = erc_params["kappa"])

trunc_pt <- 45
hist(data$fwi, freq = FALSE, breaks = 35)
curve(mix_dens_trunc_expo(x, trunc_pt, fwi_params), add = TRUE)
fwi_unif <- mix_cdf_trunc_expo(data$fwi, trunc_pt, fwi_params)

fwi_exp <- qexp(fwi_unif)
erc_exp <- qexp(erc_unif)

indices_expo_friend <- cbind(x = erc_exp, y = fwi_exp) |>
  as_tibble() |>
  mutate(
    r = x + y,
    w1 = x / r,
    w2 = y / r,
    dataset = "friendmtn"
  )

fire_data_list <- grab_top_n(indices_expo, n0 = 185, N = nrow(indices_expo))
n <- nrow(fire_data_list$cloud_tib)
fire_data_list$cloud_tib |>
  ggplot(aes(x / log(n), y / log(n), color = high)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(x = y_lb / log(n), y = x_lb / log(n), color = "red"), linewidth = 0.5) +
  theme_classic() +
  scale_color_manual(values = cols) +
  # scale_x_continuous(expand = c(0,0.01)) +
  # scale_y_continuous(expand = c(0,0.01)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.025))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  theme(
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.position = "none"
  ) +
  xlab(expression("ERC" / "log(n)")) +
  ylab(expression("FWI" / "log(n)"))

fire_data_list$cloud_tib |>
  ggplot(aes(w1, r / log(n), color = high)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(x = w1, y = r0_w / log(n), color = "red"), linewidth = 0.75) +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  theme(
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.position = "none"
  ) +
  xlab(expression("W")) +
  ylab(expression("R" / "log(n)"))

qs::qsave(fire_data_list, file = "data/raw/friendmtn_expo.qs")

## plots of transformation
erc_tibble <- tibble(idx_val = data$erc, idx = "erc") |>
  mutate(
    grid_val = seq(min(idx_val), max(idx_val), length.out = n()),
    dens_val = g1_pdf(grid_val, sigma = erc_params["sigma"], xi = erc_params["xi"], kappa = erc_params["kappa"])
  )

fwi_tibble <- tibble(idx_val = data$fwi, idx = "fwi") |>
  mutate(
    grid_val = seq(min(idx_val), max(idx_val) + 1, length.out = n()),
    dens_val = mix_dens_trunc_expo(grid_val, trunc_pt, fwi_params)
  )

friend_tibble <- rbind(erc_tibble, fwi_tibble)

erc_plot <- erc_tibble |> ggplot(aes(x = idx_val)) +
  geom_histogram(aes(y = after_stat(density)),
    bins = 25,
    boundary = 0, alpha = 0.5, position = "identity", color = "grey34", fill = "grey"
  ) +
  geom_line(aes(x = grid_val, y = dens_val), color = "blue", linewidth = 1) +
  xlab("ERC") +
  ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "none",
    legend.background = element_rect(fill = "transparent", color = "transparent")
  )

fwi_plot <- fwi_tibble |> ggplot(aes(x = idx_val)) +
  geom_histogram(aes(y = after_stat(density)),
    bins = 22,
    boundary = 0, alpha = 0.5, position = "identity", color = "grey34", fill = "grey"
  ) +
  geom_line(aes(x = grid_val, y = dens_val), color = "blue", linewidth = 1) +
  xlab("FWI") +
  ylab("Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "none",
    legend.background = element_rect(fill = "transparent", color = "transparent")
  )
erc_plot + fwi_plot

friend_transforms <- erc_plot + plot_spacer() + fwi_plot + plot_layout(widths = c(0.95, 0.05, 0.95))
friend_transforms / plot_spacer() / redstone_transforms + plot_layout(heights = c(0.95, 0.05, 0.95))

ggsave("figures/fire_hists_transform.pdf",
  dpi = 320,
  bg = "transparent",
  width = 8,
  height = 8
)
knitr::plot_crop("figures/fire_hists_transform.pdf")

# ggsave("figures/fire_hists_friendmtn.pdf",
#        dpi = 320,
#        bg = "transparent",
#        width = 8,
#        height = 4)
# knitr::plot_crop("figures/fire_hists_friendmtn.pdf")

friend_qq <- indices_expo |>
  select(x, y) |>
  rename(ERC = x, FWI = y) |>
  pivot_longer(cols = 1:2, names_to = "idx", values_to = "value") |>
  ggplot(aes(sample = value)) +
  stat_qq_band(distribution = "exp", dparams = list(rate = 1), bandType = "ks", fill = "transparent", color = "darkgrey", linewidth = 0.75) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, col = "red", linetype = "dashed") +
  stat_qq_point(distribution = "exp", dparams = list(rate = 1), alpha = 0.75) +
  scale_x_continuous(limits = c(0, 14), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  facet_grid(. ~ idx) +
  theme_classic() +
  theme(
    panel.spacing.x = unit(1, "cm"),
    panel.spacing.y = unit(0.5, "cm"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text.x = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2))
  )

friend_qq / plot_spacer() / redstone_qq + plot_layout(heights = c(0.95, 0.05, 0.95))
ggsave("figures/qq_plot_marg_transform.pdf",
  dpi = 320,
  bg = "transparent",
  width = 8,
  height = 8
)
knitr::plot_crop("figures/qq_plot_marg_transform.pdf")

indices_expo <- rbind(indices_expo_friend, indices_expo_red) |>
  select(x, y, dataset) |>
  rename(ERC = x, FWI = y) |>
  pivot_longer(cols = 1:2, names_to = "idx", values_to = "value") |>
  group_by(idx, dataset) |>
  arrange(value, .by_group = TRUE) |>
  mutate(
    empir_prob = qexp(row_number() / (n() + 1)),
    ulb = qexp(qbeta(0.025, row_number(), n() + 1 - row_number())),
    uub = qexp(qbeta(0.975, row_number(), n() + 1 - row_number()))
  ) |>
  ungroup()

indices_expo |>
  ggplot(aes(x = empir_prob, y = value)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, col = "red", linetype = "dashed") +
  geom_line(aes(x = empir_prob, y = ulb), linetype = "dashed", color = "darkgrey") +
  geom_line(aes(x = empir_prob, y = uub), linetype = "dashed", color = "darkgrey") +
  scale_x_continuous(limits = c(NA, 12), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(NA, 12), expand = expansion(mult = c(0, 0.01))) +
  facet_grid(dataset ~ idx, scales = "free") +
  theme_classic() +
  theme(
    panel.spacing.x = unit(1, "cm"),
    panel.spacing.y = unit(0.5, "cm"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text.x = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2))
  )

indices_expo |>
  ggplot(aes(sample = value)) +
  stat_qq_band(distribution = "exp", dparams = list(rate = 1), bandType = "ell", fill = "transparent", color = "darkgrey", linewidth = 0.75) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, col = "red", linetype = "dashed") +
  stat_qq_point(distribution = "exp", dparams = list(rate = 1), alpha = 0.75) +
  scale_x_continuous(limits = c(NA, 13), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(limits = c(NA, 13), expand = expansion(mult = c(0, 0.01))) +
  facet_grid(dataset ~ idx, axes = "all", axis.labels = "margins") +
  labs(
    x = "Theoretical exponential quantiles",
    y = "Sample quantiles"
  ) +
  theme_classic() +
  theme(
    panel.spacing.x = unit(0.75, "cm"),
    panel.spacing.y = unit(0.75, "cm"),
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text.x = element_text(size = rel(1.5)),
    strip.text.y = element_blank(),
    axis.text = element_text(size = rel(1.3)),
    axis.title = element_text(size = rel(1.3))
  )

ggsave("figures/qq_plot_marg_transform.pdf",
  dpi = 320,
  bg = "transparent",
  width = 10,
  height = 10
)
knitr::plot_crop("figures/qq_plot_marg_transform.pdf")


## marg transform plots all together
both_tibs <- rbind(red_tibble |> mutate(dataset = "redstone"), friend_tibble |> mutate(dataset = "friendmtn"))

both_tibs |> ggplot(aes(x = idx_val)) +
  geom_histogram(data = ~ subset(., dataset == "friendmtn" & idx == "erc"), aes(y = after_stat(density)), bins = 25, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_histogram(data = ~ subset(., dataset == "friendmtn" & idx == "fwi"), aes(y = after_stat(density)), bins = 22, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_histogram(data = ~ subset(., dataset == "redstone" & idx == "erc"), aes(y = after_stat(density)), bins = 25, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_histogram(data = ~ subset(., dataset == "redstone" & idx == "fwi"), aes(y = after_stat(density)), bins = 30, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_line(aes(x = grid_val, y = dens_val), color = "blue", linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(dataset ~ idx, scales = "free", axes = "all", axis.labels = "margins") +
  theme_classic() +
  ylab("Density") +
  xlab("Fire index") +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.2)),
    axis.title = element_text(size = rel(1.2)),
    legend.text = element_text(size = rel(1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "none",
    legend.background = element_rect(fill = "transparent", color = "transparent")
  )

red_tibble |> ggplot(aes(x = idx_val)) +
  geom_histogram(data = ~ subset(., idx == "erc"), aes(y = after_stat(density)), bins = 25, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_histogram(data = ~ subset(., idx == "fwi"), aes(y = after_stat(density)), bins = 30, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_line(aes(x = grid_val, y = dens_val), color = "blue", linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  facet_wrap(. ~ idx,
    scales = "free",
    labeller = as_labeller(c(erc = "ERC", fwi = "FWI")),
    strip.position = "bottom"
  ) +
  theme_classic() +
  ylab("Density") +
  xlab(NULL) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.3)),
    axis.title = element_text(size = rel(1.3)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = rel(1.8)),
    panel.spacing.x = unit(0.75, "cm", data = NULL)
  )

ggsave("figures/marg_transform_redstone.pdf",
  bg = "transparent",
  dpi = 320,
  width = 10,
  height = 5
)
knitr::plot_crop("figures/marg_transform_redstone.pdf")

friend_plots <- friend_tibble |> ggplot(aes(x = idx_val)) +
  geom_histogram(data = ~ subset(., idx == "erc"), aes(y = after_stat(density)), bins = 25, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_histogram(data = ~ subset(., idx == "fwi"), aes(y = after_stat(density)), bins = 22, alpha = 0.5, boundary = 0, color = "gray34", fill = "grey") +
  geom_line(aes(x = grid_val, y = dens_val), color = "blue", linewidth = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) +
  facet_wrap(. ~ idx,
    scales = "free",
    strip.position = "bottom",
    labeller = as_labeller(c(erc = NULL, fwi = NULL))
  ) +
  theme_classic() +
  ylab("Density") +
  xlab(NULL) +
  theme(
    panel.background = element_rect(fill = "transparent", color = "transparent"),
    plot.background = element_rect(fill = "transparent", color = "transparent"),
    axis.text = element_text(size = rel(1.3)),
    axis.title = element_text(size = rel(1.3)),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing.x = unit(0.75, "cm", data = NULL)
  )

friend_plots / red_plots + plot_layout(heights = c(1, 1))

ggsave("figures/marg_transform.pdf",
  bg = "transparent",
  dpi = 320,
  width = 10,
  height = 10
)
knitr::plot_crop("figures/marg_transform.pdf")
