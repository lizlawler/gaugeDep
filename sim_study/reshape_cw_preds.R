# =============================================================================
# Reshapes and saves the Campbell-Wadsworth (C-W) model predictions from the
# per-file MLE output format into tidy tibbles, one per dep_type/dep_level/box
# combination. Also computes MSE tables against the true probabilities. Output
# is consumed by combine_ls_to_others.R and combine_both_ang_to_others.R.
#
# Run interactively after all C-W MLE and prediction jobs have completed.
# Inputs:    samplers/campbell_wadsworth/mle_and_preds/{dep_type}/
# Outputs:   figures/is_preds_boxplots/joint/pred_tibbles/{dep_type}_{dep_level}_{box}_cw.qs
#            figures/is_preds_boxplots/joint/mse_tables/{dep_type}_{dep_level}_{box}_cw.qs
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(evd)
library(progressr)
library(RcppSimdJson)
library(gaugeDependence)
library(qs)
library(grafify)
library(stringr)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

# Shared ground-truth probability functions
source("sim_study/true_prob_utils.R")

# combine_fits_by_level_box --------------------------------------------------
# Reads all C-W MLE/prediction files for a given dep_type and dep_level,
# filters to a specific prediction box, and stacks into a single tibble.
combine_fits_by_level_box <- function(dep_type, dep_level, box_num) {
  fit_files <- list.files(sprintf("samplers/campbell_wadsworth/mle_and_preds/%s/", dep_type),
                          pattern    = dep_level,
                          full.names = TRUE)
  map_dfr(fit_files, ~ qread(.x)$preds) |>
    filter(box == box_num) |>
    select(pred, dataset)
}

# create_predictions_tibble --------------------------------------------------
# Assembles the C-W prediction tibble for one scenario, saves it, and computes
# an MSE table against the true probability. Box dimensions:
#   b1: X1 in (10,12), X2 in (10,12)  [symmetric high]
#   b2: X1 in (10,12), X2 in (6,8)    [asymmetric]
#   b3: X1 in (10,12), X2 in (2,4)    [highly asymmetric]
create_predictions_tibble <- function(dep_type, dep_level, box_num) {
  dim1 <- c(10, 12)
  dim2 <- switch(box_num, b1 = dim1, b2 = c(6, 8), c(2, 4))

  preds_tib <- combine_fits_by_level_box(dep_type, dep_level, box_num) |>
    mutate(ang_dens = "C-W", method = NA) |>
    rename(preds = pred)
  qsave(preds_tib, sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_cw.qs",
                            dep_type, dep_level, box_num))
  print(sprintf("C-W predictions for %s %s, in box %s have been saved", dep_level, dep_type, box_num))

  true_prob <- get_true_prob(dep_type, dep_level, dim1, dim2)

  mse_table <- preds_tib |>
    mutate(truth = true_prob, diff = preds - truth) |>
    group_by(method, ang_dens) |>
    summarise(mse = mean(diff^2), .groups = "drop") |>
    mutate(log_mse = log(mse), rmse_norm = sqrt(mse) / true_prob) |>
    arrange(rmse_norm)
  qsave(mse_table, sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_cw.qs",
                            dep_type, dep_level, box_num))
}

dep_types  <- c("gauss", "logistic", "husler_reiss")
dep_levels <- c("high", "mid", "low")
boxes      <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)

with_progress({
  p <- progressor(steps = nrow(all_combos))
  apply(all_combos, 1, function(row) {
    create_predictions_tibble(dep_type  = row["dep_types"],
                              dep_level = row["dep_levels"],
                              box_num   = row["boxes"])
    p()
  })
})
