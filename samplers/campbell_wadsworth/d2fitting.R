args <- commandArgs(trailingOnly=TRUE)
dep_type <- args[1]
dep_level <- args[2]
data_start <- args[3]
data_end <- args[4]

library(geometry)
library(geometricMVE)
library(evd)
library(mvtnorm)
library(tidyr)
library(dplyr)
library(PWLExtremes)
source("samplers/campbell_wadsworth/mod_sim2d.R")

par.locs = seq(0,1,length.out=11)
gfun <- function(w,par,locs=par.locs) gfun.2d(x=w,par=par,ref.angles=locs)

find_k <- function(box, ctau_val, radial_mle) {
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  # create fake data to use in determining k value
  pseudo_pred <- expand_grid(x1_pseudo = seq(dim1[1], dim1[2], length.out=15), 
                             x2_pseudo = seq(dim2[1], dim2[2], length.out=15)) |> 
    mutate(r_pseudo = x1_pseudo + x2_pseudo,
           w1_pseudo = x1_pseudo / r_pseudo,
           w2_pseudo = x2_pseudo / r_pseudo)
  
  # determine ideal value of k using the above
  gw_pseudo <- sapply(pseudo_pred$w1_pseudo, function(x) gfun(x, radial_mle))
  poss_k <- pseudo_pred$r_pseudo * gw_pseudo / ctau_val
  return(max(round(min(poss_k), 1) - 0.1, 1))
}

make_preds <- function(n = 50000, k = 1, box, radial_mle, angular_mle, par.locs, r, w, rexc) {
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  if(k > 1) {
    sim_list <- sim.2d.joint.mod(nsim = n, k.vals = k,
                                 gfun = gfun, par = radial_mle, fW.par = angular_mle,
                                 par.locs = par.locs, r = r,w = w)[[1]]
    xstar <- sim_list$xstar
    r_over_k <- mean(sim_list$iw)
  } else {
    xstar <- sim.2d.joint.mod(nsim = n, k.vals = 1,
                              gfun = gfun, par = radial_mle, fW.par = angular_mle,
                              par.locs = par.locs, r = r,w = w)[[1]]
    r_over_k <- 1
  }
  
  prob.est <- mean(xstar[,1] > dim1[1] & xstar[,1] < dim1[2] & xstar[,2] > dim2[1] & xstar[,2] < dim2[2]) * 
    r_over_k * length(rexc)/length(r)
  return(prob.est)
}

dep_type <- "gauss"
dep_level <- "mid"
i <- 11
for(i in data_start:data_end) {
  start <- Sys.time()
  data <- RcppSimdJson::fload(sprintf("data/%s/%s_%s.json", dep_type, dep_level, i))
  
  # obtain radii and angles
  r <- data$R
  w <- data$W
  x <- cbind(r * w, r * (1-w))
  
  # estimate the threshold
  qr <- radial.quants.L1.KDE.2d(r, w, tau=0.95, bww=0.05, bwr=0.05)
  # keep the exceedances
  r0w <-qr$r0w
  wpts <- qr$wpts
  excind<- r > r0w
  rexc <- r[excind]
  wexc <- w[excind]
  r0w <- r0w[excind]
  
  # Fit the models
  model.fit.R.bounded <- fit.pwlin.2d(r = rexc, r0w = r0w, w = wexc, locs = par.locs, pen.const=1, method="BFGS", bound.fit=T)
  model.fit.W <- fit.pwlin.2d(r = rexc, r0w = r0w, w = wexc, locs = par.locs, pen.const=NULL, fW.fit=T, method="BFGS")
  
  # make predictions
  gw_fitted <- sapply(w, function(x) gfun(x, model.fit.R.bounded$mle))
  ctau_fitted <- quantile(gw_fitted * r, 0.95)
  boxes <- c("b1", "b2", "b3")
  preds_by_box <- sapply(boxes, function(x) {
    k <- find_k(box = x, ctau_val = ctau_fitted, radial_mle = model.fit.R.bounded$mle)
    pred <- make_preds(n = 50000, k = k, box = x, 
                       radial_mle = model.fit.R.bounded$mle, angular_mle = model.fit.W$fW.mle,
                       par.locs = par.locs, r = r, w = w, rexc = rexc)
    return(tibble(k = k, box = x, pred = pred))
  }, simplify = FALSE) |> bind_rows()
  end <- Sys.time()
  
  # qsave(x = results, file = sprintf("samplers/nimble/ang_mix_mcmc_fits/%s/%s_%s.qs",
  #                                   dep_type, dep_level, i))
  # print(paste0("Successfully saved MCMC stick breaking fit for dataset number: ", i))
}
