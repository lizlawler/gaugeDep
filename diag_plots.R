library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(evd)
library(progressr)
library(RcppSimdJson)
library(gaugeDependence)
library(qs)
library(grafify)
source("extraction_scripts/extract_post_params_real_data.R")
# 
options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

data_type <- "redstone"

gauge_functions <- list(
  gauss = gauss_gauge,
  inv_log = inv_log_gauge,
  rectangular = rectangular_gauge,
  logistic = logistic_gauge,
  asym_log = asym_log_gauge,
  dirichlet = dirichlet_gauge
)

# Grab gauge function by string
get_gauge_function <- function(type_str) {
  if (!type_str %in% names(gauge_functions)) {
    stop("Unknown gauge type: ", type_str)
  }
  return(gauge_functions[[type_str]])
}

ptgamma <- function(x, xmin, alpha, beta) {
  num <- pgamma(x, shape = alpha, rate = beta, lower.tail = FALSE)
  denom <- pgamma(xmin, shape = alpha, rate = beta, lower.tail = FALSE)
  return(1 - num/denom)
}

cdf_by_gauge <- function(dataset, dataname, gauge, likelihood) {
  idx <- dataset$idx
  r <- dataset$R[idx]
  w <- dataset$W[idx]
  r0w <- dataset$r0_w[idx]
  n0 <- dataset$n0
  
  # read in posterior params
  post_radial <- extract_post_params_radial(gauge, likelihood, dataname, TRUE)
  post_radial_dep <- as.numeric(post_radial[2:length(post_radial)])
  post_radial_alpha <- post_radial[["alpha"]]
  
  # calculate rate parameter
  gauge_fn <- get_gauge_function(gauge)
  gauge_vals <- gauge_fn(w, 1 - w, post_radial_dep)
  
  return(ptgamma(r, r0w, post_radial_alpha, gauge_vals) |> 
           as_tibble() |> 
           mutate(method = gauge, 
                  id = 1:n0))
}

cdf_by_lhood <- function(dataname, likelihood) {
  data <- qs::qread(sprintf("data/%s_expo.qs", dataname))
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  return(lapply(gauge_library, function(x) cdf_by_gauge(dataset = data, 
                                                        dataname = dataname,
                                                        gauge = x,  
                                                        likelihood = likelihood)) |> 
           bind_rows())
}

weighted_cdf_by_lhood <- function(likelihood, dataname) {
  cdf <- cdf_by_lhood(dataname = dataname,
                      likelihood = likelihood)
  wts_mix <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_mix.qs", dataname, likelihood))
  wtd_cdf <- suppressMessages(cdf |> left_join(wts_mix) |>
                                mutate(stacking_preds = value * stacking,
                                       pseudo_boot = pseudobma_boot * value,
                                       pseudo_noboot = pseudobma_noboot * value) |>
                                group_by(id) |>
                                summarize(stacking_predictions = sum(stacking_preds),
                                          pseudobma_boot_preds = sum(pseudo_boot),
                                          pseudobma_noboot_preds = sum(pseudo_noboot)) |>
                                ungroup()) |> 
    pivot_longer(cols = -id, names_to = "method", values_to = "cdf") |>
    mutate(method = case_when(grepl("stacking", method) ~ 'Stacking',
                              grepl("noboot", method) ~ 'Pseudo-BMA',
                              grepl("boot", method) ~ 'Pseudo-BMA+'),
           method = as.factor(method),
           likelihood = likelihood)
  return(wtd_cdf)
}

test <- weighted_cdf_by_lhood("trunc", "redstone")
trunc_test <- test |> group_by(method) |> arrange(cdf, .by_group = TRUE) |>
  mutate(theor_prob = row_number() / (n() + 1),
         ulb = qbeta(0.025, row_number(), n() + 1 - row_number()),
         uub = qbeta(0.975, row_number(), n() + 1 - row_number())) |>
  ungroup()
trunc_test |> ggplot(aes(x = theor_prob, y = cdf, group = method, color = method)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(aes(x = theor_prob, y = ulb), linetype = "dashed", color = "darkgrey") + geom_line(aes(x = theor_prob, y = uub), linetype = "dashed", color = "darkgrey") + 
  labs(title = "truncated likelihood PP plot",
       x = "Theoretical Prob",
       y = "Empirical Probability") +
  theme_classic()

trunc_test |> ggplot(aes(x = qexp(theor_prob), y = qexp(cdf), group = method, color = method)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(aes(x = qexp(theor_prob), y = qexp(ulb)), linetype = "dashed", color = "darkgrey") + geom_line(aes(x = qexp(theor_prob), y = qexp(uub)), linetype = "dashed", color = "darkgrey") + 
  labs(title = "truncated likelihood QQ plot",
       x = "Theoretical Prob",
       y = "Empirical Probability") +
  theme_classic()

cens_test <- weighted_cdf_by_lhood("cens", "redstone")
cens_test <- cens_test |> group_by(method) |> arrange(cdf, .by_group = TRUE) |>
  mutate(theor_prob = row_number() / (n() + 1),
         ulb = qbeta(0.025, row_number(), n() + 1 - row_number()),
         uub = qbeta(0.975, row_number(), n() + 1 - row_number())) |>
  ungroup()
cens_test |> ggplot(aes(x = theor_prob, y = cdf, group = method, color = method)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(aes(x = theor_prob, y = ulb), linetype = "dashed", color = "darkgrey") + 
  geom_line(aes(x = theor_prob, y = uub), linetype = "dashed", color = "darkgrey") + 
  labs(title = "censored likelihood PP plot",
       x = "Theoretical Prob",
       y = "Empirical Probability") +
  theme_classic()

cens_test |> ggplot(aes(x = qexp(theor_prob), y = qexp(cdf), group = method, color = method)) + geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(aes(x = qexp(theor_prob), y = qexp(ulb)), linetype = "dashed", color = "darkgrey") + geom_line(aes(x = qexp(theor_prob), y = qexp(uub)), linetype = "dashed", color = "darkgrey") + 
  labs(title = "censored likelihood QQ plot",
       x = "Theoretical Prob",
       y = "Empirical Probability") +
  scale_x_continuous(limits = c(0,8)) +
  theme_classic()

