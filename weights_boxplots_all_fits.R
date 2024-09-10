library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(progressr)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

dep_level_lhood_thres <- expand_grid(dep_types = c("gauss", "logistic", "husler_reiss"),
                                     levels = c("low", "mid", "high"),
                                     lhoods = c("trunc", "cens"),
                                     thres = c("marg", "ctau")) |>
  filter(!(dep_types == "husler_reiss" & thres == "ctau"))


make_wts_df <- function(weights_file) {
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  temp <- readRDS(weights_file) |>
    bind_rows() |> 
    mutate(method = rep(gauge_library, 100)) |>
    mutate(stacking = as.numeric(stacking),
           pseudobma_boot = as.numeric(pseudobma_boot),
           pseudobma_noboot = as.numeric(pseudobma_noboot)) |>
    mutate(dataset = rep(1:100, times = rep(6, 100)))
  return(temp)
}

make_wts_boxplot <- function(dep_type, dep_level, likelihood, threshold) {
  basename <- paste0(dep_type, "_", dep_level, "_", likelihood, "_", threshold)
  weights_file <- paste0("stacking_weights/", basename, "_wts.RDS")
  temp_df <- make_wts_df(weights_file) |> pivot_longer(cols = -c(dataset, method),
                                                       names_to = 'bma_method',
                                                       values_to = 'weights') |>
    mutate(bma_method = case_when(grepl("stacking", bma_method) ~ 'Stacking',
                                  grepl("noboot", bma_method) ~ 'Pseudo-BMA',
                                  grepl("boot", bma_method) ~ 'Pseudo-BMA+'),
           bma_method = as.factor(bma_method),
           method = str_to_sentence(method),
           method = as.factor(case_when(grepl("Inv", method) ~ 'Inv. log.',
                                        grepl("Asym", method) ~ 'Asym. log.',
                                        grepl("Dirichlet", method) ~ 'Dirichlet',
                                        grepl("Rectangular", method) ~ 'Rectangular',
                                        .default = method)))
  temp_plot <- temp_df |> ggplot(aes(x = method, y = weights, fill = bma_method)) + 
    geom_boxplot() + theme_classic() + 
    theme(axis.title.x = element_text(size = rel(1.0)),
          axis.text.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(1.0)),
          axis.text.y = element_text(size = rel(0.9)),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          axis.text = element_text(size = rel(1.2)),
          axis.title = element_text(size = rel(1.2))) +
    xlab("Model") + ylab("Weights") + labs(fill = "")
  saveRDS(temp_plot, file = paste0("bma_update_deck/weights_plots/", basename, ".RDS"))
}

with_progress({
  # Create a progress handler
  p <- progressor(steps = nrow(dep_level_lhood_thres))
  
  # Apply the function using apply and update the progress bar
  apply(dep_level_lhood_thres, 1, function(row) {
    p()  # Update the progress bar
    make_wts_boxplot(dep_type = row["dep_types"], 
                     dep_level = row["levels"],
                     likelihood = row["lhoods"],
                     threshold = row["thres"]
    )
  })
})


