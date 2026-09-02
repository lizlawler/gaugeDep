# =============================================================================
# Compares the theoretical conditional exceedance probability P(R > r0(W) | W=w)
# under the true bivariate model against the Gamma(2, g(W)) approximation used
# in the proposed model, via their ratio. Validates the Gamma radial assumption
# by showing the ratio stays near 1 across the angle W (produces the manuscript
# figure faceted over Gaussian and logistic dependence at three levels each).
#
# Inputs:    (generates data internally using theoretical distributions)
# Outputs:   figures/cdf_ratio_truth_v_model.png
#            viz_scripts/data_ratio_cdf.qs (cached ratio data)
#
# The true (X, Y) densities on exponential margins, and their transforms to
# (R, W), are courtesy of Jenny Wadsworth.
# =============================================================================

library(evd)
library(mvtnorm)
library(tidyverse)
library(geometricMVE)
library(grafify)

# Find the marginal threshold q leaving n0 points above it, and return the
# lower-boundary curve of the exceedance region in (r, w) form (r0w_lb, w_lb)
# along with the exceedances themselves.
grab_top_n <- function(r, w, x, y, n0 = 1, N = 10000) {
  tau <- (N-n0)/N
  q1 <- quantile(x, tau)
  q2 <- quantile(y, tau)
  q <- max(q1, q2)
  joint <- cbind(x, y) |> as_tibble()
  idx <- which(joint$x > q | joint$y > q)
  eps <- 0.001
  while (length(idx) > n0) {
    q <- q + eps
    idx <- which(joint$x > q | joint$y > q)
  }
  r0w <- ifelse(w > 0.5, q/w, q/(1-w))
  x_lb <- ifelse(w < 0.5, q, q*y/x)
  y_lb <- ifelse(w > 0.5, q, q*x/y)
  r0w_lb <- x_lb + y_lb
  w_lb <- x_lb / r0w_lb
  r_above <- r[idx]
  w_above <- w[idx]
  return(list(x_lb = x_lb, y_lb = y_lb, 
              r0w_lb = r0w_lb, w_lb = w_lb,
              w_above = w_above, r_above = r_above))
}

# Overlay the true (black) and Gamma-approximation (red) exceedance
# probabilities against W -- a quick visual check of the two curves.
plot_both <- function(data) {
  ylimits <- c(min(c(data$true_exc_prob, data$gam_exc_prob)), max(c(data$true_exc_prob, data$gam_exc_prob)))
  return(data |> ggplot(aes(x = w_top, y = true_exc_prob)) + geom_line() + theme_classic() + ylim(ylimits) +
           geom_line(aes(x = w_top, y = gam_exc_prob), col = "red"))
}

# Plot the ratio of (1 - true) / (1 - Gamma) exceedance probabilities against
# W, with the three prediction-box angle ranges shaded. Ratio near 1 means
# the Gamma approximation matches the truth.
b1_col <- get_graf_colours("contrast_blue")
b2_col <- get_graf_colours("contrast_red")
b3_col <- get_graf_colours("contrast_yellow")

plot_ratio <- function(data) {
  p <- data |> ggplot(aes(x = w_top, y = ratio)) + geom_line(linewidth = 0.75) + theme_classic() +
    scale_x_continuous(limits = c(NA, 1), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    annotate("rect", xmin = 5/11, xmax = 6/11, ymin = -Inf, ymax = Inf, fill = b1_col, color = b1_col, alpha = 0.65) +
    annotate("rect", xmin = 5/9, xmax = 6/9, ymin = -Inf, ymax = Inf, fill = b2_col, color = b2_col, alpha = 0.65) +
    annotate("rect", xmin = 5/7, xmax = 6/7, ymin = -Inf, ymax = Inf, fill = b3_col, color = b3_col, alpha = 0.65) +
    ylab("Ratio") + xlab(expression("W"["1"])) +
    theme(panel.background = element_rect(fill='transparent', color='transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'),
          axis.text = element_text(size = rel(1.3)),
          axis.title = element_text(size = rel(1.3)),
          legend.background = element_rect(fill='transparent', color='transparent'))
  return(p)
}

##################################################################################
##                                  GAUSSIAN                                    ##
##################################################################################

# The next block builds, for the Gaussian model: the joint density on
# exponential margins, its (R, W) form, and the conditional survivor
# function P(R > r | W = w) via numerical integration (infinity ~ 35).
f_gauss_exp <- function(x, y, alpha) {
  z1 <- qnorm(1 - exp(-x))
  z2 <- qnorm(1 - exp(-y))
  
  J <- ((1 - pnorm(z1)) / dnorm(z1)) * ((1 - pnorm(z2)) / dnorm(z2))
  f <- dmvnorm(x=c(z1,z2), mean = c(0,0), sigma = matrix(c(1, alpha, alpha, 1), nrow = 2)) * J
  return(f)
}

# Density of (R,W)
#==================
f_gauss_exp <- Vectorize(f_gauss_exp, vectorize.args = c("x","y"))

f_gauss_rw <- function(r, w, alpha) {
  return(exp(log(r) + log(f_gauss_exp(x = r * w, y=r * (1-w), alpha = alpha))))
}

# Density of R|W [NB numerical integral approximates infinity with 100]
#======================================================================
f_gauss_rw <- Vectorize(f_gauss_rw, vectorize.args = "r")

f_gauss_r_given_w <- function(r, w, alpha) {
  f_gauss_rw(r, w, alpha) / integrate(f_gauss_rw, w = w, lower = 0, upper = 35, alpha = alpha)$value
}

# P(R>r|W=w) [NB numerical integral approximates infinity with 100]
#==================================================================
cond_surv_function_gauss <- function(r, w, alpha) {
  integrate(f_gauss_rw, w = w, alpha = alpha, lower = r, upper = 35)$value / integrate(f_gauss_rw, w = w, lower = 0, upper = 35, alpha = alpha)$value
}

########################################################################

# Generate data
gauss_gen <- function(dep, resamp_n) {
  x <- rmvnorm(5000, mean = c(0,0), sigma = matrix(c(1, dep, dep, 1), nrow = 2))
  x1 <- qexp(pnorm(x[,1]))
  x2 <- qexp(pnorm(x[,2]))
  r <- x1 + x2
  w <- x1 / r
  
  top_pts <- grab_top_n(r, w, x1, x2, n0 = 250)
  
  # thin out number of points because the calculation will take much longer
  length(unique(round(top_pts$w_lb, 2)))
  dupe_sets <- split(seq_along(round(top_pts$w_lb, 2)), round(top_pts$w_lb, 2)) 
  dupe_sets_sub_samp <- lapply(dupe_sets, function(list_idx) {
    if(length(list_idx) > resamp_n) {
      sample(list_idx, resamp_n, replace = FALSE)
    } else {
      list_idx
    }
  })
  sub_idx <- as.integer(unlist(dupe_sets_sub_samp))
  
  w_top <- top_pts$w_lb[sub_idx]
  r0w_top <- top_pts$r0w_lb[sub_idx]
  
  true_exc_prob <- gam_exc_prob <- NULL
  for(i in seq_along(w_top)) { 
    true_exc_prob[i] <- tryCatch(cond_surv_function_gauss(r = r0w_top[i], w = w_top[i], alpha = dep), error = function(e) NA)
    gam_exc_prob[i] <- pgamma(r0w_top[i], shape = 2, rate = gauge_gaussian(c(w_top[i], 1 - w_top[i]), par = dep), lower.tail = F)
  }
  
  data_for_plot <- cbind(w_top = w_top, true_exc_prob, gam_exc_prob) |> as_tibble() |> 
    # mutate(ratio = true_exc_prob/gam_exc_prob)
    mutate(ratio = (1-true_exc_prob)/(1-gam_exc_prob))
  return(data_for_plot)
}

gauss_high_data <- gauss_gen(0.9, resamp_n = 5)
plot_both(gauss_high_data)
plot_ratio(gauss_high_data)

gauss_mid_data <- gauss_gen(0.5, resamp_n = 4)
plot_both(gauss_mid_data)
plot_ratio(gauss_mid_data)

gauss_low_data <- gauss_gen(0.1, resamp_n = 4)
plot_both(gauss_low_data)
plot_ratio(gauss_low_data)

all_gauss_data <- gauss_high_data |> mutate(dep_type = "gauss_high") |>
  rbind(gauss_mid_data |> mutate(dep_type = "gauss_mid")) |>
  rbind(gauss_low_data |> mutate(dep_type = "gauss_low"))

##################################################################################
##                                  LOGISTIC                                    ##
##################################################################################

# Density of (X,Y) in exponential margins
#========================================
f_log_exp <- function(x, y, alpha) {
  z1 <- -1/log(1-exp(-x))
  z2 <- -1/log(1-exp(-y))
  
  # For numerical purposes:
  if(!is.finite(z1)){z1 <- exp(x)}
  if(!is.finite(z2)){z2 <- exp(y)}
  
  J <- (z1^2 * (1 - exp(-1/z1)) / exp(-1/z1)) * (z2^2 * (1-exp(-1/z2)) /exp(-1/z2))
  #  print(c(z1,z2,J))
  f <- dbvevd(x = c(z1,z2), dep = alpha, mar1 = c(1,1,1)) * J
  return(f)
}

# Density of (R,W)
#==================
f_log_exp <- Vectorize(f_log_exp, vectorize.args = c("x","y"))

f_log_rw<-function(r, w, alpha) {
  return(r * f_log_exp(x = r * w, y = r * (1-w), alpha = alpha))
}

# Density of R|W [NB numerical integral approximates infinity with 100]
#======================================================================
f_log_rw <- Vectorize(f_log_rw, vectorize.args = "r")

f_log_r_given_w<-function(r, w, alpha) {
  f_log_rw(r, w, alpha) / integrate(f_log_rw, w = w, lower = 0, upper = 100, alpha = alpha)$value
}

# P(R>r|W=w) [NB numerical integral approximates infinity with 100]
#==================================================================
cond_surv_function_log <- function(r, w, alpha) {
  integrate(f_log_rw, w = w, alpha = alpha, lower = r, upper = 100)$value / integrate(f_log_rw, w=w, lower = 0, upper = 100, alpha = alpha)$value
}

########################################################################

# Data
logistic_gen <- function(dep, resamp_n) {
  x <- rbvevd(5000, dep = dep, mar1=c(1,1,1))
  x <- qexp(exp(-1/x))
  r <- x[,1] + x[,2]
  w <- x[,1] / r
  
  # thin out number of points because the calculation will take much longer
  top_pts <- grab_top_n(r, w, x[,1], x[,2], n0 = 250)
  length(unique(round(top_pts$w_lb, 2)))
  dupe_sets <- split(seq_along(round(top_pts$w_lb, 2)), round(top_pts$w_lb, 2)) 
  dupe_sets_sub_samp <- lapply(dupe_sets, function(list_idx) {
    if(length(list_idx) > resamp_n) {
      sample(list_idx, resamp_n, replace = FALSE)
    } else {
      list_idx
    }
  })
  sub_idx <- as.integer(unlist(dupe_sets_sub_samp))
  
  w_top <- top_pts$w_lb[sub_idx]
  r0w_top <- top_pts$r0w_lb[sub_idx]
  
  true_exc_prob <- gam_exc_prob <- NULL
  for(i in seq_along(w_top)) {
    true_exc_prob[i] <- tryCatch(cond_surv_function_log(r = r0w_top[i], w = w_top[i], alpha = dep), error = function(e) NA)
    gam_exc_prob[i] <- pgamma(r0w_top[i], shape = 2, rate = gauge_rvad(c(w_top[i], 1 - w_top[i]), par = dep), lower.tail = F)
  }
  
  data_for_plot <- cbind(w_top = w_top, true_exc_prob, gam_exc_prob) |> as_tibble() |> 
    # mutate(ratio = true_exc_prob/gam_exc_prob)
    mutate(ratio = (1-true_exc_prob)/(1-gam_exc_prob))
  return(data_for_plot)
}

logistic_high_data <- logistic_gen(0.1, resamp_n = 9)
plot_both(logistic_high_data)
plot_ratio(logistic_high_data)

logistic_mid_data <- logistic_gen(0.5, resamp_n = 4)
plot_both(logistic_mid_data)
plot_ratio(logistic_mid_data)

logistic_low_data <- logistic_gen(0.9, resamp_n = 4)
plot_both(logistic_low_data)
plot_ratio(logistic_low_data)

all_logistic_data <- logistic_high_data |> mutate(dep_type = "logistic_high") |>
  rbind(logistic_mid_data |> mutate(dep_type = "logistic_mid")) |>
  rbind(logistic_low_data |> mutate(dep_type = "logistic_low"))

all_data <- rbind(all_gauss_data, all_logistic_data) |>
  mutate(dep_type = factor(dep_type, 
                           levels = c("gauss_low", "gauss_mid", "gauss_high",
                                      "logistic_low", "logistic_mid", "logistic_high")))

all_data |> ggplot(aes(x = w_top, y = ratio, group = dep_type)) + 
  scale_x_continuous(limits = c(NA, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0.85, 1.15), expand = expansion(mult = c(0, 0.01))) +
  annotate("rect", xmin = 5/11, xmax = 6/11, ymin = -Inf, ymax = Inf, fill = b1_col, color = b1_col, alpha = 0.65) +
  annotate("rect", xmin = 5/9, xmax = 6/9, ymin = -Inf, ymax = Inf, fill = b2_col, color = b2_col, alpha = 0.65) +
  annotate("rect", xmin = 5/7, xmax = 6/7, ymin = -Inf, ymax = Inf, fill = b3_col, color = b3_col, alpha = 0.65) +
  facet_wrap(. ~ dep_type, axes = "all") +
  geom_hline(yintercept = 1, col = "grey34", linetype = "dashed") +
  geom_line(linewidth = 0.75) + theme_classic() +
  ylab("Ratio") + xlab(expression("W"["1"])) +
  theme(panel.background = element_rect(fill='transparent', color='transparent'),
        panel.spacing.x = unit(0.95, "cm", data = NULL),
        panel.spacing.y = unit(0.95, "cm", data = NULL),
        plot.background = element_rect(fill='transparent', color='transparent'),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.8)),
        strip.text.x = element_blank(),
        legend.background = element_rect(fill='transparent', color='transparent'))

ggsave(filename = "figures/cdf_ratio_truth_v_model.png",
       dpi = 320,
       bg = "transparent",
       width = 12,
       height = 8)
knitr::plot_crop("figures/cdf_ratio_truth_v_model.png")
qs::qsave(all_data, "viz_scripts/data_ratio_cdf.qs")
