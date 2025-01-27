gauge_gaussian3d<-function(xyz,par)
{
  x<-xyz[1];y<-xyz[2];z<-xyz[3]
  S<-matrix(c(1,par[1:2],par[1],1,par[3],par[2],par[3],1),3,3,byrow=T)
  return(t(sqrt(c(x,y,z)))%*%solve(S)%*%sqrt(c(x,y,z)))
}


gauge_rvad_full3d<-function(xyz,par)
{
  x<-xyz[1];y<-xyz[2];z<-xyz[3]
  if(length(par)==1){
  return((x+y+z)/par+ min(x,y,z)*(1-3/par))
  } else{
    return((x/par[1]+y/par[2]+z/par[3])+ min(x,y,z)*(1-1/par[1]-1/par[2]-1/par[3]))
  }
}

gauge_clayton3d<-function(xyz,par)
{
  x<-xyz[1];y<-xyz[2];z<-xyz[3]
  return(max(x,y,z)*(1+3*par)-(x+y+z)*par)
}

# gauge_rvad_pw3d
# =================

# Allows pairwise AD in a 3d vector
# If all thetas are set to 1 then all pairs are AD (put triplewise AI)
# For only (x,y) and (x,z) AD, set theta23=10e10
# For only (x,y) and (y,z) AD, set theta13=10e10
# For only (x,z) and (y,z) AD, set theta12=10e10

gauge_rvad_pw3d<-function(xyz,par,theta12=1,theta13=1,theta23=1)
{
  x<-xyz[1];y<-xyz[2];z<-xyz[3]
  #alpha12<-par[1]; alpha13<-par[2]; alpha23<-par[3]

  k<-1
  if(theta12==1){alpha12<-par[k];k<-k+1} else{alpha12<-1}
  if(theta13==1){alpha13<-par[k];k<-k+1} else{alpha13<-1}
  if(theta23==1){alpha23<-par[k];k<-k+1} else{alpha23<-1}

  bit1<-theta12*((x+y)/alpha12 +min(x,y)*(1-2/alpha12)) + min(theta13*(z/alpha13+ min(x,z)*(1-1/alpha13)),theta23*(z/alpha23+ min(y,z)*(1-1/alpha23)))
  bit2<-theta13*((x+z)/alpha13 +min(x,z)*(1-2/alpha13)) + min(theta12*(y/alpha12+ min(x,y)*(1-1/alpha12)), theta23*(y/alpha23 + min(y,z)*(1-1/alpha23)))
  bit3<-theta23*((y+z)/alpha23 +min(y,z)*(1-2/alpha23)) + min(theta12*(x/alpha12+ min(x,y)*(1-1/alpha12)), theta13*(x/alpha13 + min(x,z)*(1-1/alpha13)))

  return(min(bit1,bit2,bit3))
}

# ================

# gauge_rvad_all3d

# Allows any combination of joint extremes in a 3d vector
# If all thetas are set to 1 then there is mass on all possible subspaces of the simplex S_2
# In order to *not* have mass on a subspace, set the corresponding theta equal to 10e10 - e.g. to not place mass on {1,2,3}, set theta123=10e10

# Order of parameters (when present): par=c(alpha12,alpha13,alpha23,alpha123)
# If e.g. theta123=10e10 then alpha123 is not needed and can be omitted
xyz <- runif(3)
xyz <- xyz / sum(xyz)
x <- xyz[1]
y <- xyz[2]
z <- xyz[3]
mxy <- min(x, y)
mxz <- min(x, z)
myz <- min(y, z)
mxyz <- min(x, y, z)
alpha1 <- runif(1)
alpha2 <- runif(1)
alpha3 <- runif(1)
alpha12 <- runif(1)
alpha13 <- runif(1)
alpha23 <- runif(1)
alpha123 <- runif(1)
theta1=1
theta2=1
theta12=1
theta13=1
theta23=1
theta123=1
# V1<-c(x1[theta1==1],
#       (x1/alpha12+mxy*(1-1/alpha12))[theta12==1])
# V2<-c(y1[theta1==1],
#       (y1/alpha12+mxy*(1-1/alpha12))[theta12==1])
# V12<-((x1+y1)/alpha12 + mxy*(1-2/alpha12))[theta12==1]
V1<-c(x[theta1==1],(x/alpha12+mxy*(1-1/alpha12))[theta12==1],(x/alpha13+mxz*(1-1/alpha13))[theta13==1],(x/alpha123+mxyz*(1-1/alpha123))[theta123==1])
V2<-c(y[theta2==1],(y/alpha12+mxy*(1-1/alpha12))[theta12==1],(y/alpha23+myz*(1-1/alpha23))[theta23==1],(y/alpha123+mxyz*(1-1/alpha123))[theta123==1])
V3<-c(z[theta3==1],(z/alpha23+myz*(1-1/alpha23))[theta23==1],(z/alpha13+mxz*(1-1/alpha13))[theta13==1],(z/alpha123+mxyz*(1-1/alpha123))[theta123==1])

V12<-c(((x+y)/alpha12 + mxy*(1-2/alpha12))[theta12==1], ((x+y)/alpha123 +mxyz*(1-2/alpha123))[theta123==1])
V13<-c(((x+z)/alpha13 + mxz*(1-2/alpha13))[theta13==1], ((x+z)/alpha123 +mxyz*(1-2/alpha123))[theta123==1])
V23<-c(((y+z)/alpha23 + myz*(1-2/alpha23))[theta23==1], ((y+z)/alpha123 +mxyz*(1-2/alpha123))[theta123==1])

V123<-((x+y+z)/alpha123 + mxyz*(1-3/alpha123))[theta123==1]


term1 <- rowSums(expand.grid(V1, V2, KEEP.OUT.ATTRS = F))
sum_term1 <- rowSums(term1)
min(sum_term1)
term2 <- V12
min(min(sum_term1, term2))
gauge_rvad_all3d<-function(xyz,par,theta1=1,theta2=1,theta3=1,theta12=1,theta13=1,theta23=1,theta123=1)
{
  x<-xyz[1];y<-xyz[2];z<-xyz[3]

  mxy<-min(x,y)
  mxz<-min(x,z)
  myz<-min(y,z)
  mxyz<-min(x,y,z)

  k<-1
  if(theta12==1){alpha12<-par[k];k<-k+1} else{alpha12<-1}
  if(theta13==1){alpha13<-par[k];k<-k+1} else{alpha13<-1}
  if(theta23==1){alpha23<-par[k];k<-k+1} else{alpha23<-1}
  if(theta123==1){alpha123<-par[k];k<-k+1} else{alpha123<-1}

  V1<-c(x[theta1==1],
        (x/alpha12+mxy*(1-1/alpha12))[theta12==1],
        (x/alpha13+mxz*(1-1/alpha13))[theta13==1],(x/alpha123+mxyz*(1-1/alpha123))[theta123==1])
  V2<-c(y[theta2==1],(y/alpha12+mxy*(1-1/alpha12))[theta12==1],(y/alpha23+myz*(1-1/alpha23))[theta23==1],(y/alpha123+mxyz*(1-1/alpha123))[theta123==1])
  V3<-c(z[theta3==1],(z/alpha23+myz*(1-1/alpha23))[theta23==1],(z/alpha13+mxz*(1-1/alpha13))[theta13==1],(z/alpha123+mxyz*(1-1/alpha123))[theta123==1])

  V12<-c(((x+y)/alpha12 + mxy*(1-2/alpha12))[theta12==1], ((x+y)/alpha123 +mxyz*(1-2/alpha123))[theta123==1])
  V13<-c(((x+z)/alpha13 + mxz*(1-2/alpha13))[theta13==1], ((x+z)/alpha123 +mxyz*(1-2/alpha123))[theta123==1])
  V23<-c(((y+z)/alpha23 + myz*(1-2/alpha23))[theta23==1], ((y+z)/alpha123 +mxyz*(1-2/alpha123))[theta123==1])

  V123<-((x+y+z)/alpha123 + mxyz*(1-3/alpha123))[theta123==1]

  # suppressWarnings is there because min(NULL) generates a warning
  suppressWarnings(term1<-min(apply(expand.grid(V1,V2,V3, KEEP.OUT.ATTRS = F),1,function(x){sum(x)})))
  suppressWarnings(term2<-min(apply(expand.grid(V1,V23, KEEP.OUT.ATTRS = F),1,function(x){sum(x)})))
  suppressWarnings(term3<-min(apply(expand.grid(V2,V13, KEEP.OUT.ATTRS = F),1,function(x){sum(x)})))
  suppressWarnings(term4<-min(apply(expand.grid(V3,V12, KEEP.OUT.ATTRS = F),1,function(x){sum(x)})))
  term5<-V123

  return(min(term1,term2,term3,term4,term5))
}

# gauge_rvad_all3d<-function(xyz,par,theta1=1,theta2=1,theta3=1,theta12=1,theta13=1,theta23=1,theta123=1)
# {
#   x<-xyz[1];y<-xyz[2];z<-xyz[3]
#
#   mxy<-min(x,y)
#   mxz<-min(x,z)
#   myz<-min(y,z)
#   mxyz<-min(x,y,z)
#
#   k<-1
#   if(theta12==1){alpha12<-par[k];k<-k+1} else{alpha12<-1}
#   if(theta13==1){alpha13<-par[k];k<-k+1} else{alpha13<-1}
#   if(theta23==1){alpha23<-par[k];k<-k+1} else{alpha23<-1}
#   if(theta123==1){alpha123<-par[k];k<-k+1} else{alpha123<-1}
#
#   V1<-c(theta1*x,theta12*(x/alpha12+mxy*(1-1/alpha12)),theta13*(x/alpha13+mxz*(1-1/alpha13)),theta123*(x/alpha123+mxyz*(1-1/alpha123)))
#   V2<-c(theta2*y,theta12*(y/alpha12+mxy*(1-1/alpha12)),theta23*(y/alpha23+myz*(1-1/alpha23)),theta123*(y/alpha123+mxyz*(1-1/alpha123)))
#   V3<-c(theta3*z,theta23*(z/alpha23+myz*(1-1/alpha23)),theta13*(z/alpha13+mxz*(1-1/alpha13)),theta123*(z/alpha123+mxyz*(1-1/alpha123)))
#
#   V12<-c(theta12*((x+y)/alpha12 + mxy*(1-2/alpha12)), theta123*((x+y)/alpha123 +mxyz*(1-2/alpha123)))
#   V13<-c(theta13*((x+z)/alpha13 + mxz*(1-2/alpha13)), theta123*((x+z)/alpha123 +mxyz*(1-2/alpha123)))
#   V23<-c(theta23*((y+z)/alpha23 + myz*(1-2/alpha23)), theta123*((y+z)/alpha123 +mxyz*(1-2/alpha123)))
#
#   V123<-theta123*((x+y+z)/alpha123 + mxyz*(1-3/alpha123))
#
#   term1<-min(apply(expand.grid(V1,V2,V3),1,function(x){sum(x)}))
#   term2<-min(apply(expand.grid(V1,V23),1,function(x){sum(x)}))
#   term3<-min(apply(expand.grid(V2,V13),1,function(x){sum(x)}))
#   term4<-min(apply(expand.grid(V3,V12),1,function(x){sum(x)}))
#   term5<-V123
#
#   return(min(term1,term2,term3,term4,term5))
# }



# Functions for mixing gauge functions additively and rescaling them to represent exponential margins
#====================================================================================================

# additivegauge.scaling.3d and additivegauge.rescale.3d allow for mixing two gauges additively;
# this could be extended in future

# The output of additivegauge.scaling.3d, ms (stands for "maxscale"), is the amount to scale
# each margin by in order to get the additive gauge to be appropriately scaled for exponential margins

#====================================================================================================

# wgrid should be a grid of points (matrix with 3 columns) covering the 2-d simplex S_2

additivegauge.scaling.3d<-function(gauge1,gauge2,par1,par2,weight,wgrid)
{
  dummy<-function(xyz)
  {
    gauge1(xyz,par=par1) + weight*gauge2(xyz,par=par2)
  }

  den<-apply(wgrid,1,dummy)
  ms1<-max(wgrid[,1]/den)
  ms2<-max(wgrid[,2]/den)
  ms3<-max(wgrid[,3]/den)

  return(c(ms1,ms2,ms3))
}

create.wgrid.3d<-function(n)
{
  wseq<-seq(0,1,len=n)
  wgrid<-expand.grid(wseq,wseq)
  wgrid<-cbind(wgrid,1-wgrid[,1]-wgrid[,2])
  neg.ind<-wgrid[,3]<0
  wgrid<-wgrid[!neg.ind,]
  return(wgrid)
}

additivegauge.rescale.3d<-function(xyz,gauge1,gauge2,par1,par2,weight,ms)
{
  xyz<-c(ms[1]*xyz[1],ms[2]*xyz[2],ms[3]*xyz[3])
  return(gauge1(xyz,par=par1) + weight*gauge2(xyz,par=par2))
}

