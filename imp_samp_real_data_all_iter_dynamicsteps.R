library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(evd)
library(progressr)
library(RcppSimdJson)
library(gaugeDependence)
library(qs)
library(grafify)
source("extraction_scripts/extract_post_params_real_data.R")

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

data_type <- "redstone"

gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")

preds_by_gauge <- function(gauge, likelihood, data, box, p, compute_preds = TRUE) {
  post_radial <- extract_post_params_radial(gauge, likelihood, data, FALSE)
  n_iter <- nrow(post_radial)  # Dynamically determine row count
  
  if (!compute_preds) {
    return(n_iter * 2)  # Return steps only
  }
  
  is_samp <- gen_is_samples(box = box)
  
  mix <- map_dbl(1:n_iter, function(i) {
    preds <- is_prob_pred(imp_samples = is_samp,
                          post_radial_params = post_radial[i, ],
                          post_angular_params = extract_post_params_ang_mix(data, FALSE)[i, ],
                          box = box,
                          gauge_type = gauge, 
                          ang_dens = "mix")
    p()  # Update progress for each iteration
    preds
  })
  
  star <- map_dbl(1:n_iter, function(i) {
    preds <- is_prob_pred(imp_samples = is_samp,
                          post_radial_params = post_radial[i, ],
                          post_angular_params = extract_post_params_ang_star(gauge, data, FALSE)[i, ],
                          box = box,
                          gauge_type = gauge, 
                          ang_dens = "star")
    p()  # Update progress for each iteration
    preds
  })
  
  preds <- cbind(mix, star) |> as_tibble() |> mutate(method = gauge, iter = 1:n_iter)
  
  return(preds)
}

preds_by_lhood <- function(likelihood, data, box, p, compute_preds = TRUE) {
  results <- lapply(gauge_library, function(x) 
    preds_by_gauge(gauge = x, likelihood = likelihood, data = data, box = box, p = p, compute_preds = compute_preds)
  )
  
  if (!compute_preds) {
    return(sum(unlist(results)))
  }
  
  return(bind_rows(results))
}

weighted_preds_by_lhood <- function(likelihood, data, box, p, compute_preds = TRUE) {
  result <- preds_by_lhood(likelihood, data, box, p, compute_preds)
  
  if (!compute_preds) {
    return(result)  # Just return step count
  }
  
  preds <- result |>
    pivot_longer(cols = c(mix, star), names_to = "ang_dens", values_to = "preds")
  
  wts_star <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_star.qs", data, likelihood)) |> mutate(ang_dens = "star")
  wts_mix <- qread(sprintf("fits_and_weights/wts_joint_model/%s_%s_mix.qs", data, likelihood)) |> mutate(ang_dens = "mix")
  wts <- rbind(wts_star, wts_mix)
  
  wtd_preds <- suppressMessages(preds |> left_join(wts) |>
                                  mutate(stacking_preds = preds * stacking,
                                         pseudo_boot = pseudobma_boot * preds,
                                         pseudo_noboot = pseudobma_noboot * preds) |>
                                  group_by(iter, ang_dens) |>
                                  summarize(stacking_predictions = sum(stacking_preds),
                                            pseudobma_boot_preds = sum(pseudo_boot),
                                            pseudobma_noboot_preds = sum(pseudo_noboot)) |>
                                  ungroup())
  
  return(wtd_preds)
}

weighted_preds <- function(data, box, p, compute_preds = TRUE) {
  lhood <- c("trunc", "cens")
  
  if (!compute_preds) {
    return(sum(sapply(lhood, function(x) 
      weighted_preds_by_lhood(likelihood = x, data = data, box = box, p = NULL, compute_preds = FALSE)
    )))
  }
  
  results <- lapply(lhood, function(x) 
    weighted_preds_by_lhood(likelihood = x, data = data, box = box, p = p, compute_preds = TRUE)
  )
  
  all_wts <- bind_rows(results)
  
  qsave(all_wts, sprintf("real_data_preds/%s_%s_all_iter.qs", data, box))
  print(sprintf("Predictions for all iterations have been saved for box: %s", box))
}

with_progress({
  all_combos <- expand_grid(data_type, boxes)
  
  total_steps <- sum(apply(all_combos, 1, function(row) {
    weighted_preds(data = row["data_type"], box = row["boxes"], compute_preds = FALSE)
  }))
  
  p <- progressr::progressor(steps = total_steps)
  
  apply(all_combos, 1, function(row) {
    weighted_preds(data = row["data_type"], box = row["boxes"], p = p, compute_preds = TRUE)
  })
})
