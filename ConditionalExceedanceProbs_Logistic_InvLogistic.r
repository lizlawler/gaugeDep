# Comparison of theoretical P(R>r_0(w)|W=w) from integrating f_R|W(r|w) 
# to P(R>r_0(w)|W=w) based on a Gamma(2,g(w)) distribution

##################################################################################
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

library(evd)
library(geometricMVE)

# Data
set.seed(123)

alpha<-0.5
x<-rbvevd(5000,dep=alpha,mar1=c(1,1,1))
x<-qexp(exp(-1/x))
r<-x[,1]+x[,2]
w<-x[,1]/r


qr<-QR.2d(r=r,w=w,method="empirical", tau=0.97) 
qrs<-QR.2d(r=r,w=w,method="smooth", tau=0.97)

# Compare true exceedance probability for r0(w) with the Gamma(2,g(w)) survival function at r0(w)
true.exc.prob<-gam.exc.prob<-NULL
for(i in 1:length(qr$wpts))
{
  true.exc.prob[i]<-cond.surv.function.log(r=qr$r.tau.wpts[i],w=qr$wpts[i],alpha=alpha)
  gam.exc.prob[i]<-pgamma(qr$r.tau.wpts[i],shape=2,rate=gauge_rvad(c(qr$wpts[i],1-qr$wpts[i]),par=alpha),lower.tail = F)
}  

# Compare true P(R>r0(w)|w) with gamma-based version
plot(qr$wpts,true.exc.prob,typ="l",ylim=c(min(c(true.exc.prob,gam.exc.prob)),max(c(true.exc.prob,gam.exc.prob))))
lines(qr$wpts,gam.exc.prob,col=2)

# Ratio -- asymptotics depends on w; quite different at w=1/2
plot(qr$wpts,true.exc.prob/gam.exc.prob,typ="l")


# Same plots but for r0(w) defined by additive quantile regression
wseq<-qrs$wpts[seq(50,950,by=10)]
r0seq<-c(qrs$r.tau.wpts[seq(50,950,by=10)])
true.exc.prob<-gam.exc.prob<-NULL

for(i in 1:length(wseq))
{
  true.exc.prob[i]<-cond.surv.function.log(r=r0seq[i],w=wseq[i],alpha=alpha)
  gam.exc.prob[i]<-pgamma(r0seq[i],shape=2,rate=gauge_rvad(c(wseq[i],1-wseq[i]),par=alpha),lower.tail = F)
}  

# True exceedance and Gamma-based version
plot(wseq,true.exc.prob,typ="l",ylim=c(min(c(true.exc.prob,gam.exc.prob)),max(c(true.exc.prob,gam.exc.prob))))
lines(wseq,gam.exc.prob,col=2)

# Ratio
plot(wseq,true.exc.prob/gam.exc.prob,typ="l")

##################################################################################
# Inverted logistic

# Density of (X,Y) in exponential margins
#========================================
f.invlog.exp<-function(x,y,alpha)
{
  z1<- 1/x
  z2<- 1/y
  
  J<- (z1^2)*(z2^2)
  f<-dbvevd(x=c(z1,z2),dep=alpha,mar1 = c(1,1,1))*J
  return(f)
}


# Density of (R,W)
#=================
f.invlog.exp<-Vectorize(f.invlog.exp,vectorize.args = c("x","y"))

f.invlog.rw<-function(r,w,alpha)
{
  return(r*f.invlog.exp(x=r*w,y=r*(1-w),alpha=alpha))
}

# Density of R|W [NB numerical integral approximates infinity with 100]
#======================================================================

f.invlog.rw<-Vectorize(f.invlog.rw,vectorize.args = "r")

f.invlog.r.given.w<-function(r,w,alpha)
{
  f.invlog.rw(r,w,alpha)/integrate(f.invlog.rw,w=w,lower = 0,upper = 100,alpha=alpha)$value
}

# P(R>r|W=w) [NB numerical integral approximates infinity with 100]
#==================================================================

cond.surv.function.invlog<-function(r,w,alpha)
{
  integrate(f.invlog.rw,w=w,alpha=alpha,lower = r,upper = 100)$value/integrate(f.invlog.rw,w=w,lower = 0,upper = 100,alpha=alpha)$value
}


########################################################################

# Now find a sequence r0(w) based on quantile regression using geometricMVE

library(evd)
library(geometricMVE)

set.seed(456)
alpha<-0.6
x<-rbvevd(5000,dep=alpha,mar1=c(1,1,1))
x<-1/x
r<-x[,1]+x[,2]
w<-x[,1]/r


qr<-QR.2d(r=r,w=w,method="empirical",tau=0.98)
qrs<-QR.2d(r=r,w=w,method="smooth",tau=0.98)



true.exc.prob<-gam.exc.prob<-NULL
for(i in 1:length(qr$wpts))
{
  true.exc.prob[i]<-cond.surv.function.invlog(r=qr$r.tau.wpts[i],w=qr$wpts[i],alpha=alpha)
  gam.exc.prob[i]<-pgamma(qr$r.tau.wpts[i],shape=2,rate=gauge_invlogistic(c(qr$wpts[i],1-qr$wpts[i]),par=1/alpha),lower.tail = F)
}  

# Plot both probabilities
plot(qr$wpts,true.exc.prob,typ="l",ylim=c(min(c(true.exc.prob,gam.exc.prob)),max(c(true.exc.prob,gam.exc.prob))))
lines(qr$wpts,gam.exc.prob,col=2)

# Ratio
plot(qr$wpts,true.exc.prob/gam.exc.prob,typ="l")


wseq<-qrs$wpts[seq(50,950,by=10)]
r0seq<-c(qrs$r.tau.wpts[seq(50,950,by=10)])
true.exc.prob<-gam.exc.prob<-NULL

for(i in 1:length(wseq))
{
  true.exc.prob[i]<-cond.surv.function.invlog(r=r0seq[i],w=wseq[i],alpha=alpha)
  gam.exc.prob[i]<-pgamma(r0seq[i],shape=2,rate=gauge_invlogistic(c(wseq[i],1-wseq[i]),par=1/alpha),lower.tail = F)
}  


plot(wseq,true.exc.prob,typ="l",ylim=c(min(c(true.exc.prob,gam.exc.prob)),max(c(true.exc.prob,gam.exc.prob))))
lines(wseq,gam.exc.prob,col=2)

plot(wseq,true.exc.prob/gam.exc.prob,typ="l")
