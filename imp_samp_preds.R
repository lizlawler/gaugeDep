library(tidyverse)
library(cmdstanr)
library(posterior)
library(evd)
source("gauge_functions_wrt_w.R")
source("gauge_functions_wrt_x.R")
library(progressr) 
library(RcppSimdJson)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

trunc_gamma <- function(x, xmin, alpha, beta) {
  unnorm_pdf <- dgamma(x, shape = alpha, rate = beta)
  norm_cst <- pgamma(xmin, shape = alpha, rate = beta, lower.tail = F)
  return(unnorm_pdf / norm_cst)
}

est_volume <- function(n = 100, pars = 0.5) {
  temp <- seq(0, 1, length.out = n)
  grid <- expand.grid(temp, temp)
  gx <- gauss_gauge_wrt_x(grid[,1], grid[,2], dep_par = pars)
  return(mean(gx <= 1))
}

dens_l1_norm <- function(w1, par_val) {
  mc_vol <- est_volume(n = 100, par_val)
  gw <- gauss_gauge(w1, par_val)
  return(1 / (gw^2 * 2 * mc_vol))
}

mix_dens <- function(w, chain_of_params) {
  alphas <- as.numeric(chain_of_params[grepl("alpha_ang", names(chain_of_params))])
  betas <- as.numeric(chain_of_params[grepl("beta_ang", names(chain_of_params))])
  weights <- as.numeric(chain_of_params[grepl("weight", names(chain_of_params))])
  n <- length(weights)
  dens <- 0.0
  for(i in 1:n) {
    dens <- dens + weights[i] * dbeta(w, alphas[i], betas[i])
  }
  return(dens)
}

true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(mvtnorm::pmvnorm(lower = c(dim1_star[1],dim2_star[1]), upper = c(dim1_star[2],dim2_star[2]), corr = corr_matrix)[1])
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

gen_is_samples <- function(box = "b1", total_n = 5000) {
  # specify dimensions of box
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  # generate samples from importance distribution
  is_samp_mvn <- mvtnorm::rmvnorm(total_n, c(mean(dim1), mean(dim2)), sigma = 2 * diag(2)) |>
    as_tibble() |>
    rename(x1 = V1, x2 = V2) |>
    filter(x1 >= 0, x2 >= 0) |> # for B3, some samples may be < 0
    mutate(r = x1 + x2, w1 = x1 / r)
  
  return(is_samp_mvn)
}

# probability calculation
is_prob_pred <- function(imp_samples, post_params, box = "b1",
                         renorm_gamma = TRUE, self_norm_is = FALSE,
                         gauge, ang_dens = "star") {
  
  # grab gauge function
  gauge_fcn <- get(paste0(gauge, "_gauge"))
  
  # specify dimensions of box
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  
  # subset IS samples based on box
  idx_in_box <- which(
    with(
      imp_samples, 
      r >= (10 / w1) & r <= (12 / w1) &
        r >= (dim2[1] / (1 - w1)) & r <= (dim2[2] / (1 - w1))
    )
  )
  
  # pull posterior parameters
  post_alpha <- post_params[["alpha"]]
  post_dep_r <- post_params[[ifelse(ang_dens == "star", "dep_r", "dep")]]
  
  # compute angular density
  ang_dens <- if (ang_dens == "star") {
    dens_l1_norm(imp_samples$w1, post_params[["dep_w"]])
  } else {
    mix_dens(imp_samples$w1, post_params)
  }
  
  # estimate RW density
  gauge_vals <- gauge_fcn(imp_samples$w1, post_dep_r)
  r0w <- qgamma(0.95, shape = post_alpha, rate = gauge_vals)
  r_giv_w_dens <- if (renorm_gamma) {
    r0w <- qgamma(0.95, shape = post_alpha, rate = gauge_vals)
    trunc_gamma(imp_samples$r, r0w, alpha = post_alpha, beta = gauge_vals)
  } else {
    dgamma(imp_samples$r, shape = post_alpha, rate = gauge_vals)
  }
  
  rw_dens <- r_giv_w_dens * ang_dens
  is_dens <- mvtnorm::dmvnorm(imp_samples[, 1:2], mean = c(mean(dim1), mean(dim2)), sigma = 2 * diag(2)) * imp_samples$r
  
  # numerator of vanilla IS and also self norm IS
  numer <- sum(rw_dens[idx_in_box] / is_dens[idx_in_box])
  # account for self normalization, if necessary
  denom <- if (self_norm_is) sum(rw_dens / is_dens) else nrow(imp_samples)
  # calculate fraction for IS part of probability calculation
  is_part <- numer / denom
  
  # # numerator of vanilla IS and also self norm IS
  # numer <- rw_dens[idx_in_box] / is_dens[idx_in_box]
  # # account for self normalization, if necessary
  # denom <- if (self_norm_is) sum(rw_dens / is_dens) else nrow(imp_samples)
  # # calculate fraction for IS part of probability calculation
  # is_wts <- numer / denom
  
  # account for renormalization of radii density
  return(is_part * ifelse(renorm_gamma, 0.05, 1))
  # return(is_wts * ifelse(renorm_gamma, 0.05, 1))
}

all_is_methods <- function(imp_samples, post_params, box = "b1",
                           gauge, ang_dens = "star") {
  vanilla_renorm <- is_prob_pred(imp_samples, post_params, box, TRUE, FALSE, gauge, ang_dens)
  vanilla_unnorm <- is_prob_pred(imp_samples, post_params, box, FALSE, FALSE, gauge, ang_dens)
  selfis_renorm <- is_prob_pred(imp_samples, post_params, box, TRUE, TRUE, gauge, ang_dens)
  selfis_unnorm <- is_prob_pred(imp_samples, post_params, box, FALSE, TRUE, gauge, ang_dens)
  return(tibble(is_method = c("vanilla_renorm", "vanilla_unnorm", "selfis_renorm", "selfis_unnorm"),
                is_pred = c(vanilla_renorm, vanilla_unnorm, selfis_renorm, selfis_unnorm)))
}


data <- fload("data/gauss/high_1.json")
W <- data$W
R <- data$R
r0w <- data$r0_w_ctau
beta_mixture_files <- list.files(path = paste0("stan/radial_angular/csv_fits/gauss/"),
                                 pattern = paste0("high", "_", 1, "_\\d{1}.csv"), full.names = TRUE)
beta_mixture_params <-  as_cmdstan_fit(beta_mixture_files) |> as_draws_df() |>
  select(any_of(contains(c("weights", "alpha","beta", "dep", ".draw", ".chain")))) |>
  group_by(.chain) |>
  summarize(across(everything(), median)) |>
  select(-.draw)
beta_mixture_params_list <- split(beta_mixture_params, beta_mixture_params$.chain)

is_samples <- gen_is_samples("b1", 5000)
wts_vec_vanilla_renorm <- is_prob_pred(is_samples, beta_mixture_params_list[[1]], box = "b1", 
                                       renorm_gamma = TRUE, self_norm_is = FALSE,
                                       gauge = "gauss", ang_dens = "beta")
wts_vec_vanilla_unnorm <- is_prob_pred(is_samples, beta_mixture_params_list[[1]], box = "b1", 
                                       renorm_gamma = FALSE, self_norm_is = FALSE,
                                       gauge = "gauss", ang_dens = "beta")
wts_vec_selfis_renorm <- is_prob_pred(is_samples, beta_mixture_params_list[[1]], box = "b1", 
                                       renorm_gamma = TRUE, self_norm_is = TRUE,
                                       gauge = "gauss", ang_dens = "beta")
wts_vec_selfis_unnorm <- is_prob_pred(is_samples, beta_mixture_params_list[[1]], box = "b1", 
                                       renorm_gamma = FALSE, self_norm_is = TRUE,
                                       gauge = "gauss", ang_dens = "beta")

hist(wts_vec_vanilla_renorm, freq = FALSE, col = "red")
hist(wts_vec_vanilla_unnorm, freq = FALSE, col = "green", add = TRUE)

hist(wts_vec_selfis_renorm / 0.05, freq = FALSE, col = "purple")
hist(wts_vec_selfis_unnorm, freq = FALSE, col = "orange", add = TRUE)

posterior_predictions <- function(gauge, level, datanum, box = "b1", true_dep) {
  # read in original data
  data <- fload(json = paste0("./data/", gauge, "/", level, "_", datanum, ".json"))
  W <- data$W
  R <- data$R
  r0w <- data$r0_w_ctau
  
  # extract parameters from angular density as mixture of betas
  beta_mixture_files <- list.files(path = paste0("stan/radial_angular/csv_fits/", gauge, "/"),
                                   pattern = paste0(level, "_", datanum, "_\\d{1}.csv"), full.names = TRUE)
  beta_mixture_params <-  as_cmdstan_fit(beta_mixture_files) |> as_draws_df() |>
    select(any_of(contains(c("weights", "alpha","beta", "dep", ".draw", ".chain")))) |>
    group_by(.chain) |>
    summarize(across(everything(), mean)) |>
    select(-.draw)
  beta_mixture_params_list <- split(beta_mixture_params, beta_mixture_params$.chain)
  
  # extract params from angular density as starshaped density
  starshaped_fit <- readRDS(paste0("mcmc_samples/", gauge, "/", level, "_", datanum, ".rds"))
  starshaped_params <- lapply(starshaped_fit, function(x) {
    alpha <- median(x$alpha[10000:25000])
    dep_w <- median(x$dep_w[10000:25000])
    dep_r <- median(x$dep_r[10000:25000])
    return(list(alpha = alpha, dep_w = dep_w, dep_r = dep_r))
  }) |> bind_rows() |> mutate(chain = row_number())
  
  is_samples <- gen_is_samples(box = box, total_n = 5000)
  all_preds <- list(starshaped = apply(starshaped_params, 1, function(row) {
    all_is_methods(is_samples, post_params = row, box = box, gauge = gauge, ang_dens = "star") |>
      mutate(chain = row[["chain"]], density = "star")
  }) |> bind_rows(), 
  beta_mix = lapply(beta_mixture_params_list, function(chain) {
    all_is_methods(is_samples, post_params = chain, box = box, gauge = gauge, ang_dens = "beta") |>
      mutate(chain = chain[[".chain"]], density = "beta")
  }) |> bind_rows()) |> 
    bind_rows() |> group_by(is_method, density) |> summarise(mean_pred = mean(is_pred))
  
  dim1 <- c(10, 12)
  dim2 <- case_when(
    box == "b1" ~ dim1,
    box == "b2" ~ c(6, 8),
    TRUE ~ c(2, 4)
  )
  truth <- if(gauge == "gauss") {
    true_gauss_prob(dim1, dim2, true_dep)
  } else if(gauge == "logistic") {
    true_bvevd_prob(dim1, dim2,true_dep, "log")
  } else {
    true_bvevd_prob(dim1, dim2,true_dep, "hr")
  }
  all_preds <- all_preds |> 
    mutate("truth" = truth, 
           "dep_level" = level,
           "box" = box,
           "dataset" = datanum,
           "dep_type" = gauge)
  return(all_preds)
}

test <- posterior_predictions("gauss", "high", 1, "b1", 0.9)

all_combos <- expand_grid(dep_type = c("gauss", "logistic"), 
                          levels = c("low","mid", "high"), 
                          datanum = 1:100,
                          boxes = c("b1", "b2", "b3")) |>
  filter(!(dep_type == "gauss" & levels == "low" & datanum %in% c(32:39, 52:59, 72:79, 95:100)),
         !(dep_type == "gauss" & levels == "mid" & datanum %in% 97:100),
         !(dep_type == "gauss" & levels == "high" & datanum %in% c(35:39, 56:59, 73:79, 94:100)),
         !(dep_type == "logistic" & levels == "low" & datanum %in% c(36:39, 54:59, 74:79, 95:100))) |>
  mutate(true_val = case_when(dep_type == "gauss" & levels == "low" ~ 0.1,
                              dep_type == "gauss" & levels == "high" ~ 0.9,
                              dep_type == "logistic" & levels == "low" ~ 0.9,
                              dep_type == "logistic" & levels == "high" ~ 0.1,
                              levels == "mid" ~ 0.5))


all_preds <- with_progress({
  # Create a progress handler
  p <- progressor(steps = nrow(all_combos))
  
  # Apply the function using apply and update the progress bar
  apply(all_combos, 1, function(row) {
    p()  # Update the progress bar
    posterior_predictions(gauge = row["dep_type"], level = row["levels"],
                          datanum = as.numeric(row["datanum"]), box = row["boxes"],
                          true_dep = as.numeric(row["true_val"]))
  })
})

all_preds_tib <- all_preds |> 
  bind_rows()
all_preds_tib <- all_preds_tib |> separate_wider_delim(cols = 'is_method', delim = "_", names = c("is_type", "norm_type"))
all_preds_tib <- all_preds_tib |> mutate(is_type = as.factor(is_type))
# |> 
#   pivot_longer(cols = c(star, beta), names_to = "angular_density", values_to = "preds") |>
#   mutate(angular_density = as.factor(angular_density))

make_boxplot <- function(tibble, gauge, level, box_num) {
  plot_title <- paste0(gauge, ", ", level, ", ", box_num)
  sub_tib <- tibble |> filter(dep_type == gauge, dep_level == level, box == box_num)
  p <- sub_tib |>
    ggplot(aes(x = norm_type, y = mean_pred, fill = density)) + 
    geom_boxplot() +
    geom_hline(yintercept = unique(sub_tib$truth), col = "darkgrey", linetype = "longdash") +
    facet_wrap(. ~ is_type) +
    theme_classic() +
    ggtitle(plot_title) +
    xlab("Normalization") + ylab("Prediction probabilities") + labs(fill = "")
  ggsave(paste0("boxplots_pred_probs/", gauge, "_", level, "_", box_num, "_unnorm_gamma.pdf"),
         plot = p,
         bg = 'transparent',
         width = 8,
         height = 7,
         dpi = 320)
  return(p)
}

make_boxplot(all_preds_tib, "gauss", "high", "b1")
make_boxplot(all_preds_tib, "gauss", "high", "b2")
make_boxplot(all_preds_tib, "gauss", "high", "b3")

make_boxplot(all_preds_tib, "gauss", "mid", "b1")
make_boxplot(all_preds_tib, "gauss", "mid", "b2")
make_boxplot(all_preds_tib, "gauss", "mid", "b3")

make_boxplot(all_preds_tib, "gauss", "low", "b1")
make_boxplot(all_preds_tib, "gauss", "low", "b2")
make_boxplot(all_preds_tib, "gauss", "low", "b3")


make_boxplot(all_preds_tib, "logistic", "high", "b1")
make_boxplot(all_preds_tib, "logistic", "high", "b2")
make_boxplot(all_preds_tib, "logistic", "high", "b3")

make_boxplot(all_preds_tib, "logistic", "mid", "b1")
make_boxplot(all_preds_tib, "logistic", "mid", "b2")
make_boxplot(all_preds_tib, "logistic", "mid", "b3")

make_boxplot(all_preds_tib, "logistic", "low", "b1")
make_boxplot(all_preds_tib, "logistic", "low", "b2")
make_boxplot(all_preds_tib, "logistic", "low", "b3")


plot <- make_boxplot(all_preds_tib, "gauss", "high", "b3")
plot + coord_cartesian(ylim = c(0, 1e-6))
all_plots <- expand_grid(types = c("gauss", "logistic"), 
                         levels = c("low", "mid", "high"),
                         boxes = c("b1", "b2", "b3"))

apply(all_plots, 1, function(row) {
  make_boxplot(all_preds_tib, gauge = row["types"], level = row["levels"], box_num = row["boxes"])
})
