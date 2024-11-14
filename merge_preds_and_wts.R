library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(patchwork)
options(rlib_name_repair_verbosity = "quiet")


layout1 <- c(
  area(t = 3, b = 5, l = 3, r = 8),
  area(t = 3, b = 5, l = 9, r = 14),
  area(t = 6, b = 8, l = 3, r = 8),
  area(t = 6, b = 8, l = 9, r = 14),
  area(t = 1, b = 2, l = 1, r = 4)
)
plot(layout1)

layout2 <- c(
  area(t = 3, b = 5, l = 3, r = 8),
  area(t = 3, b = 5, l = 9, r = 14),
  area(t = 1, b = 2, l = 1, r = 4)
)
plot(layout2)

layout3 <- c(
  area(t = 3, b = 5, l = 3, r = 8),
  area(t = 3, b = 5, l = 9, r = 14),
  area(t = 1, b = 2, l = 1, r = 4),
  area(t = 1, b = 2, l = 5, r = 8)
)
plot(layout3)


read_files <- function(filepath, object_name) {
  temp <- readRDS(filepath)
  assign(object_name, temp, envir = .GlobalEnv)
  rm(temp)
  gc()
}

read_preds_boxplots <- function(dep_type, fitted = FALSE) {
  preds_boxplots <- list.files(paste0("boxplots_pred_probs/", 
                                      dep_type, "/rds_files/"), full.names = TRUE)
  preds_boxplots <- preds_boxplots[grepl("k_not1", preds_boxplots)]
  if(dep_type != "husler_reiss") {
    preds_boxplots <- preds_boxplots[grepl("true", preds_boxplots)]
    pred_names <- str_remove(basename(preds_boxplots), "_k_not1_true_threshold_preds_boxplot_with_wc.RDS")
  } else {
    pred_names <- str_remove(basename(preds_boxplots), "_k_not1_fitted_threshold_preds_boxplot_with_wc.RDS")
  }
  
  return(for(i in seq_along(pred_names)) {
    read_files(preds_boxplots[i], pred_names[i])
  })
}

read_wts_boxplots <- function(dep_type) {
  wts_boxplots <- list.files("figures/boxplots_stacking_wts/rds_files/", 
                             pattern = dep_type,
                             full.names = TRUE)
  wts_names <- str_remove(str_remove(basename(wts_boxplots), "_boxplot.rds"), paste0(dep_type, "_"))
  
  return(for(i in seq_along(wts_names)) {
    read_files(wts_boxplots[i], wts_names[i])
  })
}


viz_preds_wts <- function(dep_type, level, box, zoom = NULL, plot = FALSE, magnify = FALSE, mag_lim = NULL) {
  if(box == "b1") {
    dims <- "(10,12) x (10,12)"
  } else if(box == "b2") {
    dims <- "(10,12) x (6,8)"
  } else {
    dims <- "(10,12) x (2,4)"
  }
  
  new_title <- paste0(str_to_title(str_replace(dep_type, "_", "-")), ", ", level, ", B = ", dims)
  
  if(dep_type != "husler_reiss") {
    p <- 
      (get(paste0(level, "_cens_ctau_wts")) +
         ggtitle("Censored, ctau") +
         theme(plot.title = element_text(hjust = 1))) +
      (get(paste0(level, "_cens_marg_wts")) + 
         ggtitle("Censored, marg") + 
         theme(plot.title = element_text(hjust = 1))) + 
      (get(paste0(level, "_trunc_ctau_wts")) +
         ggtitle("Truncated, ctau") +
         theme(plot.title = element_text(hjust = 1)) ) +
      (get(paste0(level, "_trunc_marg_wts")) + 
         ggtitle("Truncated, marg") +
         theme(plot.title = element_text(hjust = 1))) + 
      (get(paste0(level, "_", box)) + 
         theme(legend.position = "none") + 
         ggtitle(new_title) +
         coord_cartesian(ylim = ylimits)) +
      plot_layout(design = layout1, guides = "collect")
  } else if(dep_type == "husler_reiss" & magnify) {
    p <- 
      (get(paste0(level, "_cens_marg_wts")) + 
         ggtitle("Censored, marg") + 
         theme(plot.title = element_text(hjust = 1))) + 
      (get(paste0(level, "_trunc_marg_wts")) + 
         ggtitle("Truncated, marg") +
         theme(plot.title = element_text(hjust = 1))) + 
      (get(paste0(level, "_", box)) + 
         theme(legend.position = "none") + 
         ggtitle(new_title) +
         coord_cartesian(ylim = zoom)) +
      (get(paste0(level, "_", box)) + 
         theme(legend.position = "none") + 
         ggtitle("magnified to show non-zero truth") +
         coord_cartesian(ylim = mag_lim)) +
      plot_layout(design = layout3, guides = "collect")
  } else {
    p <- 
      (get(paste0(level, "_cens_marg_wts")) + 
         ggtitle("Censored, marg") + 
         theme(plot.title = element_text(hjust = 1))) + 
      (get(paste0(level, "_trunc_marg_wts")) + 
         ggtitle("Truncated, marg") +
         theme(plot.title = element_text(hjust = 1))) + 
      (get(paste0(level, "_", box)) + 
         theme(legend.position = "none") + 
         ggtitle(new_title) +
         coord_cartesian(ylim = zoom)) +
      plot_layout(design = layout2, guides = "collect")
  }
  
  ggsave(filename = paste0("bma_update_deck/preds_and_wts/", dep_type, "_", level, "_", box, "_preds_wts.pdf"), 
         p,
         height = 8,
         width = 14,
         dpi = 320)
  if(plot) {
    return(p)
  }
}

read_preds_boxplots("husler_reiss", fitted = TRUE)
read_wts_boxplots("husler_reiss")

viz_preds_wts("husler_reiss", "low", "b1", zoom = c(0, 2.5e-8), plot=TRUE)

combos <- expand_grid(levels = c("low", "mid", "high"),
                      boxes = c("b1", "b2", "b3"))
apply(combos, 1, function(row) viz_preds_wts("husler_reiss", row["levels"], row["boxes"]))

