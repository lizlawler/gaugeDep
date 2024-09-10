library(cmdstanr)
library(MCMCvis)
library(posterior)
library(tidyverse)
library(patchwork)

extract_posterior_est <- function(gauge, dep_level, data_num) {
  csvfiles <- list.files(path = paste0("stan/angular/csv_fits/", gauge, "/"), 
                         pattern = paste0(dep_level, "_", data_num, "_v2_\\d{1}.csv"),
                         full.names = TRUE)
  fit <- as_cmdstan_fit(csvfiles) |> as_draws_df() |>
    select(-c(".iteration", ".chain")) |>
    rename(draw = ".draw") |> 
    select(dep)
  # if (gauge != "dirichlet") {
  #   return(fit |> as_draws_df() |>
  #            select(-c(".iteration", ".chain")) |>
  #            rename(draw = ".draw") |> 
  #            select(dep))
  # } else {
  #   return(fit |> as_draws_df() |>
  #            select(-c(".iteration", ".chain")) |>
  #            rename(draw = ".draw") |> 
  #            select(theta1, theta2))
  # }
  return(median(fit$dep))
}


plot_boxplot <- function(gauge, dep_level, true_dep) {
  temp <- sapply(1:100, function(x) extract_posterior_est(gauge, dep_level, x))
  temp |> as_tibble() |>
    ggplot(aes(y = value, x = 0)) + 
    geom_boxplot(fill = "grey", width = 0.5) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      axis.text = element_text(size = rel(1.2)),
      axis.title = element_text(size = rel(1.2))) +
    geom_hline(yintercept = true_dep, col = "blue", linetype = 2, linewidth = 1.2) +
    ylim(0, 1) +
    xlim(-0.5, 0.5)
}

patch_boxplots <- function(gauge) {
  if(gauge == "gauss") {
    levels <- list(low = 0.1, mid = 0.5, high = 0.9)
  } else {
    levels <- list(low = 0.9, mid = 0.5, high = 0.1)
  }
  high_dep_plot <- plot_boxplot(gauge, "high", true_dep = levels$high)
  mid_dep_plot <- plot_boxplot(gauge, "mid", true_dep = levels$mid)
  low_dep_plot <- plot_boxplot(gauge, "low", true_dep = levels$low)
  
  patched_plot <- high_dep_plot | mid_dep_plot | low_dep_plot
  return(patched_plot)
}

gauss_plots <- patch_boxplots("gauss")
logistic_plots <- patch_boxplots("logistic")

ggsave(filename = "figures/angular_dep_boxplot_logistic.pdf",
       logistic_plots,
       bg = "transparent",
       width = 15,
       height = 10,
       dpi = 320)
