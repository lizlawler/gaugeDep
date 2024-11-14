library(evd)
library(mvtnorm)
library(tidyverse)
library(geometricMVE)

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


# Comparison of theoretical P(R>r_0(w)|W=w) from integrating f_R|W(r|w) 
# to P(R>r_0(w)|W=w) based on a Gamma(2,g(w)) distribution

##################################################################################
# Gaussian

# Density of (X,Y) in exponential margins
#========================================
f.gauss.exp<-function(x,y,alpha)
{
  
  z1 <- qnorm(1 - exp(-x))
  z2 <- qnorm(1 - exp(-y))

  J <- ((1 - pnorm(z1)) / dnorm(z1)) * ((1 - pnorm(z2)) / dnorm(z2))
  f <- dmvnorm(x=c(z1,z2), mean = c(0,0), sigma = matrix(c(1, alpha, alpha, 1), nrow = 2))*J
  return(f)
}

# Density of (R,W)
#==================
f.gauss.exp<-Vectorize(f.gauss.exp,vectorize.args = c("x","y"))

f.gauss.rw<-function(r,w,alpha)
{
  # top <- w + (1 - w) - 2 * alpha * sqrt(w * (1 - w))
  # gw <- top / (1 - alpha^2)
  # alpha_exp <- gw / 2
  return(exp(log(r) + log(f.gauss.exp(x=r*w,y=r*(1-w),alpha=alpha))))
}

# Density of R|W [NB numerical integral approximates infinity with 100]
#======================================================================

f.gauss.rw<-Vectorize(f.gauss.rw,vectorize.args = "r")

f.gauss.r.given.w<-function(r,w,alpha)
{
  f.gauss.rw(r,w,alpha)/integrate(f.gauss.rw,w=w,lower = 0, upper = 35, alpha=alpha)$value
}

# P(R>r|W=w) [NB numerical integral approximates infinity with 100]
#==================================================================

cond.surv.function.gauss<-function(r,w,alpha)
{
  integrate(f.gauss.rw,w=w,alpha=alpha,lower = r,upper = 35)$value/integrate(f.gauss.rw,w=w,lower = 0, upper = 35, alpha=alpha)$value
}

########################################################################

# Now find a sequence r0(w) based on quantile regression using geometricMVE

# Data
set.seed(123)

alpha<-0.9
x<- rmvnorm(5000, mean = c(0,0), sigma = matrix(c(1, alpha, alpha, 1), nrow = 2))
x1 <- qexp(pnorm(x[,1]))
x2 <- qexp(pnorm(x[,2]))
r<-x1+x2
w<-x1/r

top_pts <- grab_top_n(r, w, x1, x2, n0 = 250)
length(unique(round(top_pts$w_lb, 2)))
dupe_sets <- split(seq_along(round(top_pts$w_lb, 2)), round(top_pts$w_lb, 2))
dupe_sets_sub_samp <- lapply(dupe_sets, function(list_idx) {
  if(length(list_idx) > 4) {
    sample(list_idx, 4, replace = FALSE)
  } else {
    list_idx
  }
})
sub_idx <- as.integer(unlist(dupe_sets_sub_samp))
plot(w, r, pch = 20)
points(top_pts$w_lb[sub_idx], top_pts$r0w_lb[sub_idx], pch = 20, col = "blue")
w_top <- top_pts$w_lb[sub_idx]
r0w_top <- top_pts$r0w_lb[sub_idx]

true.exc.prob<-gam.exc.prob<-NULL
for(i in seq_along(w_top)) { ## exclude w values (and associated thresholds) that may cause the integral to diverge
  true.exc.prob[i] <- tryCatch(cond.surv.function.gauss(r = r0w_top[i], w = w_top[i], alpha=alpha), error=function(e) NA)
  gam.exc.prob[i] <- pgamma(r0w_top[i], shape = 2, rate = gauge_gaussian(c(w_top[i],1-w_top[i]), par = alpha), lower.tail = F)
}

# Compare true P(R>r0(w)|w) with gamma-based version
data_for_plot <- cbind(w = w_top, true.exc.prob, gam.exc.prob) |> as_tibble() |> mutate(ratio = true.exc.prob/gam.exc.prob)

ylimits <- c(min(c(true.exc.prob,gam.exc.prob)),max(c(true.exc.prob,gam.exc.prob)))
data_for_plot |> ggplot(aes(x = w_top, y = true.exc.prob)) + geom_line() + theme_classic() + ylim(ylimits) +
  geom_line(aes(x = w_top, y = gam.exc.prob), col = "red")

# Ratio -- asymptotics depends on w; quite different at w=1/2
data_for_plot |> ggplot(aes(x = w_top, y = ratio)) + geom_line() + theme_classic() +
  geom_vline(xintercept = 5/9, col = "red") +
  geom_vline(xintercept = 6/9, col = "red") + 
  geom_vline(xintercept = 5/11, col = "blue") +
  geom_vline(xintercept = 6/11, col = "blue") +
  geom_vline(xintercept = 5/7, col = "green") +
  geom_vline(xintercept = 6/7, col = "green") 


####################################################################################
# Logistic

# Density of (X,Y) in exponential margins
#========================================
f.log.exp<-function(x,y,alpha)
{
  z1<- -1/log(1-exp(-x))
  z2<- -1/log(1-exp(-y))
  
  # For numerical purposes:
  if(!is.finite(z1)){z1<-exp(x)}
  if(!is.finite(z2)){z2<-exp(y)}
  
  J<- (z1^2*(1-exp(-1/z1))/exp(-1/z1))*(z2^2*(1-exp(-1/z2))/exp(-1/z2))
  #  print(c(z1,z2,J))
  f<-dbvevd(x=c(z1,z2),dep=alpha,mar1 = c(1,1,1))*J
  return(f)
}

# Density of (R,W)
#==================
f.log.exp<-Vectorize(f.log.exp,vectorize.args = c("x","y"))

f.log.rw<-function(r,w,alpha)
{
  return(r*f.log.exp(x=r*w,y=r*(1-w),alpha=alpha))
}

# Density of R|W [NB numerical integral approximates infinity with 100]
#======================================================================

f.log.rw<-Vectorize(f.log.rw,vectorize.args = "r")

f.log.r.given.w<-function(r,w,alpha)
{
  f.log.rw(r,w,alpha)/integrate(f.log.rw,w=w,lower = 0,upper = 100,alpha=alpha)$value
}

# P(R>r|W=w) [NB numerical integral approximates infinity with 100]
#==================================================================

cond.surv.function.log<-function(r,w,alpha)
{
  integrate(f.log.rw,w=w,alpha=alpha,lower = r,upper = 100)$value/integrate(f.log.rw,w=w,lower = 0,upper = 100,alpha=alpha)$value
}

########################################################################

# Now find a sequence r0(w) based on quantile regression using geometricMVE

# Data
set.seed(123)

alpha<-0.9
x<-rbvevd(5000,dep=alpha,mar1=c(1,1,1))
x<-qexp(exp(-1/x))
r<-x[,1]+x[,2]
w<-x[,1]/r

top_pts <- grab_top_n(r, w, x[,1], x[,2], n0 = 250)
length(unique(round(top_pts$w_lb, 2)))
dupe_sets <- split(seq_along(round(top_pts$w_lb, 2)), round(top_pts$w_lb, 2))
dupe_sets_sub_samp <- lapply(dupe_sets, function(list_idx) {
  if(length(list_idx) > 3) {
    sample(list_idx, 3, replace = FALSE)
  } else {
    list_idx
  }
})
sub_idx <- as.integer(unlist(dupe_sets_sub_samp))
plot(w, r, pch = 20)
points(top_pts$w_lb[sub_idx], top_pts$r0w_lb[sub_idx], pch = 20, col = "blue")
w_top <- top_pts$w_lb[sub_idx]
r0w_top <- top_pts$r0w_lb[sub_idx]

true.exc.prob<-gam.exc.prob<-NULL
for(i in seq_along(w_top)) { ## exclude w values (and associated thresholds) that may cause the integral to diverge
  true.exc.prob[i] <- tryCatch(cond.surv.function.log(r = r0w_top[i], w = w_top[i], alpha=alpha), error=function(e) NA)
  gam.exc.prob[i] <- pgamma(r0w_top[i], shape = 2, rate = gauge_rvad(c(w_top[i],1-w_top[i]), par = alpha), lower.tail = F)
}

# Compare true P(R>r0(w)|w) with gamma-based version
data_for_plot <- cbind(w = w_top, true.exc.prob, gam.exc.prob) |> as_tibble() |> mutate(ratio = true.exc.prob/gam.exc.prob)

ylimits <- c(min(c(true.exc.prob,gam.exc.prob)),max(c(true.exc.prob,gam.exc.prob)))
data_for_plot |> ggplot(aes(x = w_top, y = true.exc.prob)) + geom_line() + theme_classic() + ylim(ylimits) +
  geom_line(aes(x = w_top, y = gam.exc.prob), col = "red")

# Ratio -- asymptotics depends on w; quite different at w=1/2
data_for_plot |> ggplot(aes(x = w_top, y = ratio)) + geom_line() + theme_classic() +
  geom_vline(xintercept = 5/9, col = "red") +
  geom_vline(xintercept = 6/9, col = "red") + 
  geom_vline(xintercept = 5/11, col = "blue") +
  geom_vline(xintercept = 6/11, col = "blue") +
  geom_vline(xintercept = 5/7, col = "green") +
  geom_vline(xintercept = 6/7, col = "green") 
