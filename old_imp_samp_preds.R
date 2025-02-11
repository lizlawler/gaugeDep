

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

make_boxplot <- function(tibble, gauge, level, box) {
  plot_title <- paste0(gauge, ", ", level, ", ", box)
  sub_tib <- tibble |> filter(dep_type == gauge, dep_level == level, box == box)
  p <- sub_tib |>
    ggplot(aes(x = norm_type, y = mean_pred, fill = density)) + 
    geom_boxplot() +
    geom_hline(yintercept = unique(sub_tib$truth), col = "darkgrey", linetype = "longdash") +
    facet_wrap(. ~ is_type) +
    theme_classic() +
    ggtitle(plot_title) +
    xlab("Normalization") + ylab("Prediction probabilities") + labs(fill = "")
  ggsave(paste0("boxplots_pred_probs/", gauge, "_", level, "_", box, "_unnorm_gamma.pdf"),
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
  make_boxplot(all_preds_tib, gauge = row["types"], level = row["levels"], box = row["boxes"])
})
