library(geometricMVE)
library(evd)

data <- rbvevd(500,dep = 0.1,model = "hr")
plot(data,pch=20)
x <- qexp(pgev(data[,1], loc = 0, scale = 1, shape = 0))
y <- qexp(pgev(data[,2], loc = 0, scale = 1, shape = 0))
plot(x/log(500),y/log(500),pch=20)


# Example use:

# d=2

# set.seed(1)
# n<-2500
# x<-rbvevd(n, dep=0.5, mar1=c(0,1,0))
# x<-qexp(exp(-exp(-x)))

# Read in data you used to create BMA predictions
gauss_high1 <- RcppSimdJson::fload("data/gauss/high_1.json")

# Radial-angular variables

# r<-x[,1]+x[,2]
# w<-x[,1]/r
r <- gauss_high1$R
w <- gauss_high1$W
x <- r*w
y <- r * (1-w)
gauss_gauge <- function(w, dep_par = 0.5) {
  top <- w + (1 - w) - 2 * dep_par * sqrt(w * (1 - w))
  return(top/(1-dep_par^2))
}
gw <- gauss_gauge(w, 0.9)
plot(w/gw, (1-w)/gw,col="red",pch=20)
# Establishing a high threshold for R|W

# Empirical:
qr<-QR.2d(r=r,w=w, method = "empirical")

# Smooth:
qrs<-QR.2d(r=r,w=w, method = "smooth")

# visualization of empirial gauge function from this
plotempiricalgauge.2d(wpts=qr$wpts,r.tau.wpts = qr$r.tau.wpts, add = T)
plotempiricalgauge.2d(wpts=qrs$wpts,r.tau.wpts = qrs$r.tau.wpts, add=T)
points(w/gw, (1-w)/gw,col="red",pch=20)

# Fitting truncated gamma model, using R exceedances defined through rolling-windows quantiles

?fit.geometricMVE.2d

excind<-r>qr$r0w
rexc<-r[excind]
wexc<-w[excind]

# Threshold value corresponding to each w:

r0w<-qr$r0w[excind]

# Check for and remove any NAs - small number can occur because of empirical gauge estimation
na.ind<-which(is.na(rexc))
na.ind

# Run the following if na.ind has any entries:
if(length(na.ind)>0){
  rexc<-rexc[-na.ind]
  wexc<-wexc[-na.ind]
  r0w<-r0w[-na.ind]}

# Fit using different basic parametric gauge functions

fit1<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,init.val=c(2,0.5),gfun=gauge_rvad) # logistic-type gauge
fit2<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,init.val=c(2,0.5),gfun=gauge_gaussian) # gaussian-type gauge
fit3<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,init.val=c(2,0.5),gfun=gauge_invlogistic) # inv-logistic type gauge
fit4<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,init.val=c(2,0.5),gfun=gauge_square) # alternative square-type gauge

fit1
fit2
fit3
fit4

min(fit1$nllh, fit2$nllh, fit3$nllh, fit4$nllh)

plotempiricalgauge.2d(wpts=qr$wpts,r.tau.wpts = qr$r.tau.wpts)
plotfittedgauge.2d(wpts=seq(0,1,len=1000),gfun=gauge_rvad, par=fit1$mle[2],add=T)
plotfittedgauge.2d(wpts=seq(0,1,len=1000),gfun=gauge_gaussian, par=fit2$mle[2],add=T)
plotfittedgauge.2d(wpts=seq(0,1,len=1000),gfun=gauge_invlogistic, par=fit3$mle[2],add=T)
plotfittedgauge.2d(wpts=seq(0,1,len=1000),gfun=gauge_square, par=fit4$mle[2],add=T)


# Fit using additively mixed gauge function

fitmix<-fit.geometricMVE.2d(r=rexc,w=wexc,r0w=r0w,init.val=c(2,0.5,0.5,1),
                            add.gauge = TRUE, gauge1 = gauge_rvad, gauge2 = gauge_gaussian,
                            par1ind = 2, par2ind = 3)


plotfittedgauge.2d(wpts=seq(0,1,len=1000),add.gauge=T,gauge1 = gauge_rvad, gauge2 = gauge_gaussian,
                   par1ind = 2, par2ind = 3, par=fitmix$mle[2:4], add=T, col=4)


# Diagnostic plots for the best fit (lowest AIC)

ppdiag.2d(r=rexc,w=wexc,r0w=r0w,gfun=gauge_square, par=fit4$mle)
qqdiag.2d(r=rexc,w=wexc,r0w=r0w,gfun=gauge_square, par=fit4$mle)

# Simulation of new pseudo-observations using the model structure

newx<-sim.2d(w=wexc,r0w=r0w,gfun=gauge_square, par=fit4$mle, nsim=10000)

# Simulation above a higher threshold
newx2<-sim.2d(w=wexc,r0w=r0w,k=2,gfun=gauge_square, par=fit4$mle, nsim=10000)

plot(x, y,pch=20)
points(newx,pch=20,col=4)
points(newx2,pch=20,col=3)

#=====================================================

# d=3

# Multivariate Gaussian data
library(mvtnorm)
set.seed(1)
n<-5000
x<-rmvnorm(n,sigma=matrix(c(1,0.3,0.5,0.3,1,0.7,0.5,0.7,1),3,3,byrow = T))
x<-qexp(pnorm(x))

r<-apply(x,1,sum)
w<-x/r

qr<-QR.3d(r = r, w=w, tau=0.95)
plotempiricalgauge.3d(qr$wpts,qr$r.tau.wpts)


excind<-r>qr$r0w
rexc<-r[excind]
wexc<-w[excind,]
# Threshold value corresponding to each w:

r0w<-qr$r0w[excind]


na.ind<-which(is.na(rexc))
na.ind

if(length(na.ind)>0){
  rexc<-rexc[-na.ind]
  wexc<-wexc[-na.ind,]
  r0w<-r0w[-na.ind]}


# Fit using different basic parametric gauge functions

fit1<-fit.geometricMVE.3d(r = rexc, w=wexc, r0w = r0w, gfun = gauge_rvad_full3d,
                          init.val = c(3,0.5),
                          lower.limit = c(0,0), upper.limit = c(100,1),
                          control=list(reltol=1e-6))

fit2<-fit.geometricMVE.3d(r = rexc, w=wexc, r0w = r0w, gfun = gauge_gaussian3d,
                          init.val = c(3,rep(0.5,len=3)),
                          control=list(reltol=1e-6))

fit3<-fit.geometricMVE.3d(r = rexc, w=wexc, r0w = r0w, gfun = gauge_clayton3d,
                          init.val = c(3,0.5),
                          control=list(reltol=1e-6))

fit1
fit2
fit3

# Compare empirical and fitted gauges
plotempiricalgauge.3d(qr$wpts,qr$r.tau.wpts)
plotfittedgauge.3d(grid = 30,gfun = gauge_gaussian3d,par = fit2$mle[2:4],add = T)

# PP and QQ plots
ppdiag.3d(r = rexc,w = wexc,r0w = r0w,par = fit2$mle,gfun = gauge_gaussian3d)
qqdiag.3d(r = rexc,w = wexc,r0w = r0w,par = fit2$mle,gfun = gauge_gaussian3d)


# Simulate and plot new data above two thresholds

xstar<-sim.3d(w = wexc,r0w = r0w,nsim = 10000,par = fit2$mle, gfun = gauge_gaussian3d)
xstar2<-sim.3d(w = wexc,r0w = r0w,nsim = 10000,par = fit2$mle, gfun = gauge_gaussian3d, k=2)

library(rgl)
plot3d(x)
points3d(xstar,col=2)
points3d(xstar2,col=3)