library(evd)
library(mvtnorm)
library(tidyverse)
library(gaugeDependence)

gauge_functions <- list(
  gauss = gauss_gauge,
  inv_log = inv_log_gauge,
  rectangular = rectangular_gauge,
  logistic = logistic_gauge,
  asym_log = asym_log_gauge,
  dirichlet = dirichlet_gauge
)

# Function to get a gauge function by name
get_gauge_function <- function(type_str) {
  if (!type_str %in% names(gauge_functions)) {
    stop("Unknown gauge type: ", type_str)
  }
  return(gauge_functions[[type_str]])
}

polar_euc_tib <- function(x, y, gauge_type, dep) {
  euc_tib <- cbind(x,y) |> as_tibble()
  gauge_fn <- get_gauge_function(gauge_type)
  polar_euc_gw <- euc_tib |> mutate(r = x+y,
                                    w = x/r,
                                    gw = pmap_dbl(list(x = w, y = 1-w, dep = dep), 
                                                  function(x, y, dep) gauge_fn(x, y, dep)))
  # gw_r = pmap_dbl(list(x = gw, y = r), function(x, y) x*y))
  return(polar_euc_gw)
}

# for Husler-Reiss data only (doesn't use gauge function threshold)
polar_euc_tib_hr <- function(x, y) {
  euc_tib <- cbind(x,y) |> as_tibble()
  polar_euc <- euc_tib |> mutate(r = x+y,
                                 w = x/r)
  return(polar_euc)
}

## 0.25, 2, 6

## Copula functions for varying dependence structures -----------
gauss <- function(N = 5000, dep = 0.5) {
  x <- rmvnorm(N, mean = c(0,0), sigma = matrix(c(1, dep, dep, 1), nrow = 2))
  u1 <- pnorm(x[,1])
  u2 <- pnorm(x[,2])
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, "gauss", dep))
}

logistic <- function(N = 5000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, "logistic", dep))
}

inv_log <- function(N = 5000, dep = 0.5) {
  x <- rbvevd(N, dep = dep, model = "log", mar1=c(1,1,1))
  y <- 1/x
  x <- y[,1]
  y <- y[,2]
  return(polar_euc_tib(x, y, "inv_log", dep))
}

asym_log <- function(N = 5000, dep = 0.5, t1 = 0.5, t2 = 0.5) {
  x <- rbvevd(N, dep = dep, model = "alog", asy = c(t1, t2))
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, "asym_log", dep))
}

dirichlet <- function(N = 5000, theta1, theta2 = 2) {
  x <- rbvevd(N, alpha = theta1, beta = theta2, model = 'ct')
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib(x, y, "dirichlet", dep))
}

husler_reiss <- function(N = 5000, dep = 1) {
  x <- rbvevd(N, dep = dep, model = "hr")
  u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
  x <- qexp(u1)
  y <- qexp(u2)
  return(polar_euc_tib_hr(x, y))
}


# husler_reiss <- function(N = 5000, dep = 1) {
#   x <- rbvevd(N, dep = dep, model = "hr")
#   u1 <- pgev(x[,1], loc = 0, scale = 1, shape = 0)
#   u2 <- pgev(x[,2], loc = 0, scale = 1, shape = 0)
#   x <- qexp(u1)
#   y <- qexp(u2)
#   return(cbind(x, y) |> as_tibble())
# }
# # 
# low_hr <- husler_reiss(dep = 0.25) |> mutate(dep = "dep = 0.25")
# mid_hr <- husler_reiss(dep = 2) |> mutate(dep = "dep = 2")
# high_hr <- husler_reiss(dep = 6) |> mutate(dep = "dep = 6")
# all_hr <- rbind(low_hr, mid_hr, high_hr) |> mutate(dep = as.factor(dep))

# ggplot(all_hr, aes(x=x, y=y,  color = dep)) + 
#   geom_point() + 
#   facet_wrap(. ~ dep, axes = "all", axis.labels = "all_x") +
#   theme_classic() +
#   theme(legend.position = "none",
#         panel.background = element_rect(fill='transparent'),
#         plot.background = element_rect(fill='transparent', color='transparent')) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0))

# 
# ggsave("bma_update_deck/hr_set.pdf",
#        height = 3.5,
#        width = 10.5,
#        bg = 'transparent',
#        dpi = 320)
# 
# plot(high_hr, pch = 20)
# points(mid_hr, pch = 20, col = "blue")
# points(low_hr, pch = 20, col = "orange")

## Creation of stan data list -------------
grab_top_n <- function(cloud_tib, n0 = 1, N = 5000) {
  tau <- (N-n0)/N
  q1 <- quantile(cloud_tib$x, tau)
  q2 <- quantile(cloud_tib$y, tau)
  q <- max(q1, q2)
  idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  eps <- 0.001
  while (length(idx) > n0) {
    q <- q + eps
    idx <- which(cloud_tib$x > q | cloud_tib$y > q)
  }
  r0_w <- cloud_tib |> mutate(r0_w = ifelse(w > 0.5, q/w, q/(1-w))) |> select(r0_w)
  if("gw_r" %in% colnames(cloud_tib)) {
    ctau <- as.numeric(quantile(cloud_tib$gw_r, tau))
    r0_w_ctau <- ctau / cloud_tib$gw
    idx_ctau <- which(cloud_tib$r > r0_w_ctau)
    return(list(q = q,
                idx = idx,
                n0 = length(idx),
                N = N,
                R = cloud_tib$r,
                W = cloud_tib$w,
                r0_w = r0_w$r0_w,
                ctau = ctau,
                r0_w_ctau = r0_w_ctau,
                n0_ctau = length(idx_ctau),
                idx_ctau = idx_ctau))
  } else {
    return(list(q = q,
                idx = idx,
                n0 = length(idx),
                N = N,
                R = cloud_tib$r,
                W = cloud_tib$w,
                r0_w = r0_w$r0_w))
  }
}

data_list <- function(N = 5000, n0 = 1, dep, cop_func) {
  og_data <- cop_func(N, dep)
  return(grab_top_n(og_data, n0 = n0, N = N))
}

gen_data_file <- function(N = 5000, n0 = 1, dep, cop_func, dep_level, iter) {
  temp <- data_list(N, n0, dep, cop_func)
  cmdstanr::write_stan_json(temp, sprintf("data/%s/%s_%s.json",
                                          deparse(substitute(cop_func)), dep_level, iter))
  # cmdstanr::write_stan_json(temp, paste0("data/", dep_type, "/", dep_level, "_", iter, ".json"))
}

dep_levels <- list(c(0.25, "low"), c(2, "mid"), c(6, "high"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:200) {
    gen_data_file(5000, 250, dep, husler_reiss, level, i)
  }
}

dep_levels <- list(c(0.1, "low"), c(0.5, "mid"), c(0.9, "high"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:200) {
    gen_data_file(5000, 250, dep, gauss, level, i)
  }
}

dep_levels <- list(c(0.1, "high"), c(0.5, "mid"), c(0.9, "low"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:200) {
    gen_data_file(5000, 250, dep, logistic, level, i)
  }
}

dep_levels <- list(c(3, "high"), c(1, "mid"), c(0.5, "low"))
for (j in seq_along(dep_levels)) {
  dep <- as.numeric(dep_levels[[j]][1])
  level <- dep_levels[[j]][2]
  for ( i in 1:200) {
    gen_data_file(5000, 1000, dep, dirichlet, level, i)
  }
}
