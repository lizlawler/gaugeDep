library(qs)
library(RcppSimdJson)
library(extRemes)
library(scales)

#### GAUSS ####

# high
data <- fload("data/gauss/high_1.json")
w <- data$W
idx <- data$idx
hist(w[idx], freq = FALSE, xlim = c(0,1))
hist(w, freq = FALSE, xlim = c(0,1), add = TRUE, col = alpha("blue", 0.5))
     qqplot(w[idx], w)

# mid
data <- fload("data/gauss/mid_1.json")
w <- data$W
idx <- data$idx
hist(w[idx], freq = FALSE, xlim = c(0,1), breaks = 30)
hist(w, freq = FALSE, xlim = c(0,1), breaks = 30, add = TRUE, col = alpha("blue", 0.5))
qqplot(w[idx], w)

# low
data <- fload("data/gauss/low_1.json")
w <- data$W
idx <- data$idx
hist(w[idx], freq = FALSE, xlim = c(0,1), breaks = 30)
hist(w, freq = FALSE, xlim = c(0,1), breaks = 30, add = TRUE, col = alpha("blue", 0.5))
qqplot(w[idx], w)

#### LOGISTIC ####

# high
data <- fload("data/logistic/high_1.json")
w <- data$W
idx <- data$idx
hist(w[idx], freq = FALSE, xlim = c(0,1))
hist(w, freq = FALSE, xlim = c(0,1), add = TRUE, col = alpha("blue", 0.5))
qqplot(w[idx], w)

# mid
data <- fload("data/logistic/mid_1.json")
w <- data$W
idx <- data$idx
hist(w[idx], freq = FALSE, xlim = c(0,1), breaks =25)
hist(w, freq = FALSE, xlim = c(0,1), breaks = 25, add = TRUE, col = alpha("blue", 0.5))
qqplot(w[idx], w)

# low
data <- fload("data/logistic/low_1.json")
w <- data$W
idx <- data$idx
hist(w[idx], freq = FALSE, xlim = c(0,1), breaks = 30)
hist(w, freq = FALSE, xlim = c(0,1), breaks = 30, add = TRUE, col = alpha("blue", 0.5))
qqplot(w[idx], w)
