library(dirichletprocess)
library(tidyverse)
library(RcppSimdJson)

## gaussian fits ----------
# low dependence
low_1_data <- fload("data/angular/gauss/low_1.json")
w1 <- low_1_data$w1
hist(w1, freq = FALSE, breaks = 30)
dp_gauss_low <- DirichletProcessBeta(w1, 1)
dp_gauss_low <- Fit(dp_gauss_low, 2000)

# mid dependence
mid_1_data <- fload("data/angular/gauss/mid_1.json")
w1 <- mid_1_data$w1
hist(w1, freq = FALSE, breaks = 30)
dp_gauss_mid <- DirichletProcessBeta(w1, 1)
dp_gauss_mid <- Fit(dp_gauss_mid, 2000)

# high dependence
high_1_data <- fload("data/angular/gauss/high_1.json")
w1 <- high_1_data$w1
hist(w1, freq = FALSE, breaks = 30)
dp_gauss_high <- DirichletProcessBeta(w1, 1)
dp_gauss_high <- Fit(dp_gauss_high, 2000)

## logistic fits ----------
# low dependence
low_1_data <- fload("data/angular/logistic/low_1.json")
w1 <- low_1_data$w1
hist(w1, freq = FALSE, breaks = 30)
dp_logistic_low <- DirichletProcessBeta(w1, 1)
dp_logistic_low <- Fit(dp_logistic_low, 2000)

# mid dependence
mid_1_data <- fload("data/angular/logistic/mid_1.json")
w1 <- mid_1_data$w1
hist(w1, freq = FALSE, breaks = 30)
dp_logistic_mid <- DirichletProcessBeta(w1, 1)
dp_logistic_mid <- Fit(dp_logistic_mid, 2000)

# high dependence
high_1_data <- fload("data/angular/logistic/high_1.json")
w1 <- high_1_data$w1
hist(w1, freq = FALSE, breaks = 30)
dp_logistic_high <- DirichletProcessBeta(w1, 1)
dp_logistic_high <- Fit(dp_logistic_high, 2000)



xGrid <- seq(0, 1, by=0.001)
postSamples <- data.frame(replicate(100, PosteriorFunction(dp_beta)(xGrid)))
postFrame <- data.frame(x=xGrid, y=rowMeans(postSamples))
lines(x = postFrame$x, y = postFrame$y)

plot(dp_beta)

# ClusterTraceplot(dp_beta)

params_list <- dp_gauss_high$clusterParameters

alpha_post <- rep(NA, length(params_list[[1]]))
beta_post <- rep(NA, length(params_list[[2]]))
for(i in seq_along(params_list[[1]])) {
  alpha_post[i] <- params_list[[1]][i] * params_list[[2]][i]
  beta_post[i] <- (1-params_list[[1]][i]) * params_list[[2]][i]
}

curve(dbeta(x, alpha_post[1], beta_post[1]), ylim = c(0,10))
for(i in 2:6) {
  curve(dbeta(x, alpha_post[i], beta_post[i]), add = TRUE)
}


bernstein_dens <- function(w, weights) {
  k <- length(weights)
  pdf <- 0.0
  for(j in 1:k) {
    pdf = pdf + weights[j] * dbeta(w, shape1 = j, shape2 = (k - j + 1))
  }
  return(pdf)
}
