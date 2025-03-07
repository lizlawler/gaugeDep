library(tidyverse)
library(qs)
library(evd)

true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(mvtnorm::pmvnorm(lower = c(dim1_star[1],dim2_star[1]), 
                          upper = c(dim1_star[2],dim2_star[2]), 
                          corr = corr_matrix)[1])
}

true_bvevd_prob <- function(dim1, dim2, dep, model_type) {
  dim1_star <- qgev(pexp(dim1), loc = 0, scale = 1, shape = 0)
  dim2_star <- qgev(pexp(dim2), loc = 0, scale = 1, shape = 0)
  upper_right <- pbvevd(q = c(dim1_star[2], dim2_star[2]), model = model_type, dep = dep)
  upper_left <- pbvevd(q = c(dim1_star[1], dim2_star[2]), model = model_type, dep = dep)
  lower_right <- pbvevd(q = c(dim1_star[2], dim2_star[1]), model = model_type, dep = dep)
  lower_left <- pbvevd(q = c(dim1_star[1], dim2_star[1]), model = model_type, dep = dep)
  return(upper_right - upper_left - lower_right + lower_left)
}

true_gauss_prob(c(10,12), c(2,4), 0.9)
levels_list <- list(low = 0.25, mid = 2, high = 6)
true_bvevd_prob(dim1 = c(10,12), dim2 = c(2,4), dep = 0.25, model_type = "hr")

plot_wc <- qread("figures/is_preds_boxplots/joint/plot_objects/logistic_low_wc_b2.qs")
plot_wc + coord_cartesian(ylim = c(0, 2.5e-5))

plot_og <- qread("figures/is_preds_boxplots/joint/plot_objects/husler_reiss_low_b3.qs")
plot_og + coord_cartesian(ylim = c(0, 1e-5))


true_bvevd_prob(c(10,12), c(2,4), 0.5, "log")
plot <- qread("figures/is_preds_boxplots/joint/plot_objects/logistic_mid_b3.qs")
plot + coord_cartesian(ylim = c(0, 2e-7))


plot <- qread("figures/is_preds_boxplots/joint/plot_objects/logistic_mid_b3.qs")
plot + coord_cartesian(ylim = c(0, 1.5e-6))
