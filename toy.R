library(evd)
library(mvtnorm)
library(tidyverse)

gauss_copula <- function(n = 1000, rho = 0.5, qmarg1, qmarg2) {
  xy <- rmvnorm(n, mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  u1 <- pnorm(xy[,1])
  u2 <- pnorm(xy[,2])
  u1u2 <- cbind(u1, u2)
  v <- qmarg1(u1)
  w <- qmarg2(u2)
  vw <- cbind(v, w)
  return(list("xy" = xy, "u1u2" = u1u2, "vw" = vw))
}
q1 <- function(x) qexp(x)
q2 <- function(x) qexp(x)

n <- 100000
low_corr_exp <- gauss_copula(n, 0.1,q1, q2)
mid_corr_exp <- gauss_copula(n, 0.5,q1, q2)
high_corr_exp <- gauss_copula(n, 0.9,q1, q2)

p1 <- low_corr_exp$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.1", " expo margins")))
ggExtra::ggMarginal(p1, type = "histogram")

p2 <- mid_corr_exp$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.5", " expo margins")))
ggExtra::ggMarginal(p2, type = "histogram")

p3 <- high_corr_exp$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.9", " expo margins")))
ggExtra::ggMarginal(p3, type = "histogram")

q1 <- function(x) qunif(x)
q2 <- function(x) qunif(x)

low_corr_unif <- gauss_copula(n, 0.1,q1, q2)
mid_corr_unif <- gauss_copula(n, 0.5,q1, q2)
high_corr_unif <- gauss_copula(n, 0.9,q1, q2)

p4 <- low_corr_unif$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.1,", " unif margins")))
ggExtra::ggMarginal(p4, type = "histogram")

p5 <- mid_corr_unif$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.5", " unif margins")))
ggExtra::ggMarginal(p5, type = "histogram")

p6 <- high_corr_unif$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.9", " unif margins")))
ggExtra::ggMarginal(p6, type = "histogram")

q1 <- function(x) qlogis(x)
q2 <- function(x) qlogis(x)

low_corr_logis <- gauss_copula(n, 0.1,q1, q2)
mid_corr_logis <- gauss_copula(n, 0.5,q1, q2)
high_corr_logis <- gauss_copula(n, 0.9,q1, q2)

p7 <- low_corr_logis$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.1,", " logis margins")))
ggExtra::ggMarginal(p7, type = "histogram")

p8 <- mid_corr_logis$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.5", " logis margins")))
ggExtra::ggMarginal(p8, type = "histogram")

p9 <- high_corr_logis$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste(rho, "=0.9", " logis margins")))
ggExtra::ggMarginal(p9, type = "histogram")


logistic_copula <- function(n = 1000, r = 0.5, qmarg1, qmarg2) {
  xy <- rbvevd(n, dep = r, model = "log")
  u1 <- pgev(xy[,1], loc = 0, scale = 1, shape = 0)
  u2 <- pgev(xy[,2], loc = 0, scale = 1, shape = 0)
  u1u2 <- cbind(u1, u2)
  v <- qmarg1(u1)
  w <- qmarg2(u2)
  vw <- cbind(v, w)
  return(list("xy" = xy, "u1u2" = u1u2, "vw" = vw))
}

n <- 100000
high_dep_exp <- logistic_copula(n, 0.1,q1, q2)
mid_dep_exp <- logistic_copula(n, 0.5,q1, q2)
low_dep_exp <- logistic_copula(n, 0.9,q1, q2)

p <- high_dep_exp$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste("r=0.1", " expo margins")))
ggExtra::ggMarginal(p, type = "histogram")

p <- mid_dep_exp$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste("r=0.5", " expo margins")))
ggExtra::ggMarginal(p, type = "histogram")

p <- low_dep_exp$vw %>% as_tibble() %>% mutate(v = v/log(n), w = w/log(n)) %>%
  ggplot(aes(v, w)) + geom_point(alpha=0.5, color = "blue") + theme_classic() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  xlab("V") + ylab("W") + ggtitle(expression(paste("r=0.9", " expo margins")))
ggExtra::ggMarginal(p, type = "histogram")


n <- 1000
dep.par <- 0.5
x <- evd::rbvevd(n, dep=dep.par, mar1=c(1,1,1))
y <- 1/x
plot(y/log(n))


