library(cmdstanr)
library(posterior)
library(tidyverse)
library(RcppSimdJson)
library(mvtnorm)
library(evd)
library(qs)
options(rlib_name_repair_verbosity = "quiet")

# create function to reshape previously extracted stacking weights -------
# make_wts_df <- function(weights_file) {
#   gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
#   temp <- readRDS(weights_file) |>
#     bind_rows() |> 
#     mutate(gauge = rep(gauge_library, 100)) |>
#     mutate(stacking = as.numeric(stacking),
#            pseudobma_boot = as.numeric(pseudobma_boot),
#            pseudobma_noboot = as.numeric(pseudobma_noboot)) |>
#     mutate(dataset = rep(1:100, times = rep(6, 100))) 
#   return(temp)
# }

make_wts_boxplot <- function(weights_file) {
  temp <- qread(weights_file) |> rename(gauge = method) |>
    pivot_longer(cols = -c("gauge", "dataset"), names_to = "method", values_to = "weights") |>
    mutate(method = case_when(grepl("stacking", method) ~ 'Stacking',
                              grepl("noboot", method) ~ 'Pseudo-BMA',
                              grepl("boot", method) ~ 'Pseudo-BMA+'),
           method = as.factor(method),
           gauge = case_when(gauge == "logistic" ~ "Logistic",
                             gauge == "gauss" ~ "Gaussian",
                             gauge == "inv_log" ~ "Inv. log.",
                             gauge == "asym_log" ~ "Asym. log.",
                             gauge == "dirichlet" ~ "Dirichlet",
                             gauge == "rectangular" ~ "Rectangular"))
  plot_name <- stringr::str_to_sentence(str_replace_all(str_remove(basename(weights_file), ".qs") , "_", ", "))
  p <- ggplot(temp, aes(x = gauge, y = weights, fill = method)) + geom_boxplot() +
    theme_classic() +
    ggtitle(plot_name) +
    xlab("Gauge function") + ylab("Stacking weights") + labs(fill = "BMA method")
  pdf_name <- paste0("figures/boxplots_stacking_wts/",
                     str_remove(basename(weights_file), ".qs"),
                     "_boxplot.pdf")
  
  qs_name <- paste0("figures/boxplots_stacking_wts/qs_files/",
                     str_remove(basename(weights_file), ".qs"),
                     "_boxplot.qs")
  ggsave(pdf_name,
         plot = p,
         dpi = 320,
         bg = "transparent",
         width = 15, height = 8)
  print(paste0(pdf_name, " has been written to disk."))
  qsave(p, file = qs_name)
  print(paste0(qs_name, " has been written to disk."))
}

make_wts_boxplot("fits_and_weights/wts_joint_model/gauss_cens_sb_high.qs")

files_to_loop <- list.files(path = "fits_and_weights/wts_joint_model/", pattern = ".qs", full.names = TRUE)
# files_to_loop <- files_to_loop[!grepl("_wc_", files_to_loop)]
sapply(files_to_loop, function(x) make_wts_boxplot(x))
