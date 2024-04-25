library(bench)
library(RcppSimdJson)
library(ggplot2)
library(doParallel)
library(cmdstanr)

Rcpp::sourceCpp("gauge_mcmc.cpp")

a <- 4
b <- 2

w <- seq(0, 1, length.out = 500)
gw <- gauss_gauge(w, 0.5)
r <- 2.59

gamma_boost_seed <- function(alpha, beta, seed) {
  set.seed(seed)
  xgamma <- rgamma(1, alpha, rate = beta)
  return(gamma_boost(xgamm, alpha, beta))
}

gamma_R_seed <- function(alpha, beta, seed) {
  set.seed(seed)
  xgamma <- rgamma(1, alpha, rate = beta)
  return(gamma_R(xgamma, alpha, beta))
}

gamma_lccdf_boost_vals <- function(alpha, gw_val, r) {
  val <- rep(NA, length(gw_val))
  for(i in seq_along(gw_val)) {
    val[i] <- gamma_lccdf_boost(alpha, gw_val[i], r)
  }
  return(val)
}

gamma_lccdf_R_vals <- function(alpha, gw_val, r) {
  val <- rep(NA, length(gw_val))
  for(i in seq_along(gw_val)) {
    val[i] <- gamma_lccdf_R(alpha, gw_val[i], r)
  }
  return(val)
}

gamma_R_seed <- function(alpha, beta, seed) {
  set.seed(seed)
  xgamma <- rgamma(1, alpha, rate = beta)
  return(gamma_R(alpha, beta, xgamma))
}

bm <- mark(
    iterations = 100,
    sample_boost = gamma_lccdf_boost_vals(a, gw, r),
    sample_R = gamma_lccdf_R_vals(a, gw, r)
)

bm |> 
  group_by(expression = expression %>% as.character()) |>
  summarise(mean_total = mean(total_time), 
            mean_mem = mean(mem_alloc),
            median_time = median(median))

autoplot(bm)
boost_v_R <- benchmark("boost" = {
  gamma_boost(4, 2, rgamma(1, 4,rate=2))
})


model <- cmdstan_model("stan/bivar_trunc_gauss.stan", compile = TRUE)

liz_v_stan_trunc_low <- benchmark("stan" = {
  model$sample(data = "data/gauss/low_10.json", 
               chains = 3,
               parallel_chains = 3)
},
"liz" = {
  idx <- fload("data/gauss/low_10.json")$idx
  W <- fload("data/gauss/low_10.json")$W
  R <- fload("data/gauss/low_10.json")$R
  r0_w <- fload("data/gauss/low_10.json")$r0_w
  set.seed(4, kind = "L'Ecuyer-CMRG") 
  registerDoParallel(cores = 4)
  trunc_chains <- foreach(t = 1:3) %dopar% {
    trunc_mcmc_mh(n_iter = 10000, W = W[idx], R = R[idx], r0_w[idx], step_size = c(0.1, 0.2), update_int = 20)
  }
  stopImplicitCluster()
},
replications = 10,
columns = c("test", "replications", "elapsed",
            "relative"))

model <- cmdstan_model("stan/bivar_cens_gauss.stan", compile = TRUE)

liz_v_stan_cens_low <- benchmark("stan" = {
  model$sample(data = "data/gauss/low_10.json", 
               chains = 3,
               parallel_chains = 3)
},
"liz" = {
  idx <- fload("data/gauss/low_10.json")$idx
  W <- fload("data/gauss/low_10.json")$W
  R <- fload("data/gauss/low_10.json")$R
  r0_w <- fload("data/gauss/low_10.json")$r0_w
  set.seed(4, kind = "L'Ecuyer-CMRG") 
  registerDoParallel(cores = 4)
  cens_chains <- foreach(t = 1:3) %dopar% {
    cens_mcmc_mh(n_iter = 5000, W, R, r0_w, step_size = c(0.05, 0.05), update_int = 20)
  }
  stopImplicitCluster()
},
replications = 10,
columns = c("test", "replications", "elapsed",
            "relative"))
