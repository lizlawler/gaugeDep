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
library(stringr)

options(rlib_name_repair_verbosity = "quiet")
handlers("cli")

true_gauss_prob <- function(dim1, dim2, dep) {
  dim1_star <- qnorm(pexp(dim1))
  dim2_star <- qnorm(pexp(dim2))
  corr_matrix <- matrix(c(1, dep, dep, 1), nrow = 2)
  return(mvtnorm::pmvnorm(lower = c(dim1_star[1],dim2_star[1]), 
                          upper = c(dim1_star[2],dim2_star[2]), 
                          corr = corr_matrix)[1])
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

create_combined_preds <- function(dep_type, dep_level, box) {
  
  dim1 <- c(10, 12)
  if(box == "b1") {
    dim2 <- dim1
  } else if(box == "b2") {
    dim2 <- c(6, 8)
  } else {
    dim2 <- c(2, 4)
  }
  
  ls_preds <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s.qs", dep_type, dep_level, box))
  bezels_preds <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_bezels.qs", dep_type, dep_level, box))
  full_preds <- rbind(ls_preds, bezels_preds)
  
  # determine true probability
  if(dep_type == "gauss") {
    levels_list <- list(low = 0.1, mid = 0.5, high = 0.9, high_wc = 0.8)
    true_prob <- true_gauss_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]))
  } else if(dep_type == "logistic"){
    levels_list <- list(low = 0.9, mid = 0.5, high = 0.1, low_wc = 0.8, mid_wc = 0.4)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "log")
  } else {
    levels_list <- list(low = 0.1, mid = 1, high = 3)
    true_prob <- true_bvevd_prob(dim1 = dim1, dim2 = dim2, dep = as.numeric(levels_list[dep_level]), model_type = "hr")
  }
  
  # create boxplot
  plot <- full_preds |> ggplot(aes(x = ang_dens, y = pmax(preds, .Machine$double.eps), fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = pmax(true_prob, .Machine$double.eps), col = "darkgrey", linetype = "longdash") +
    scale_y_log10() +
    ggtitle(sprintf("%s, %s, (%s) x (%s)", dep_type, dep_level, paste(dim1, collapse = ","), paste(dim2, collapse = ","))) + 
    theme_classic() +
    theme(panel.background = element_rect(fill='transparent', color='transparent'),
          plot.background = element_rect(fill='transparent', color='transparent'),
          axis.text = element_text(size = rel(1.2)),
          axis.title = element_text(size = rel(1.2)),
          legend.text = element_text(size = rel(1.2)),
          legend.title = element_text(size = rel(1.2)),
          legend.background = element_rect(fill='transparent', color='transparent')) +
    scale_fill_grafify(breaks = ~ .x[!is.na(.x)], palette = "r4", ColSeq = FALSE) + # use colorblind-friendly palette
    labs(fill = "BMA method") +
    scale_x_discrete(labels=c("BezELS", "Cens., Mix", "Cens., Star", "Trunc., Mix", "Trunc., Star")) +
    xlab("Likelihood, Angular Density") + ylab("Prediction probabilities")
  ggsave(sprintf("figures/is_preds_boxplots/joint/%s_%s_%s_withBez.pdf", dep_type, dep_level, box),
         plot = plot,
         bg = 'transparent',
         width = 8,
         height = 7,
         dpi = 320)
  qsave(plot, sprintf("figures/is_preds_boxplots/joint/plot_objects/%s_%s_%s_withBez.qs",dep_type, dep_level, box))
  print(sprintf("Joint boxplot for %s %s, in box %s has been saved", dep_level, dep_type, box))
  
  # combine MSE tables
  ls_mse <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s.qs", dep_type, dep_level, box))
  bez_mse <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_bezels.qs", dep_type, dep_level, box))
  full_mse <- rbind(ls_mse, bez_mse) |> arrange(mse)
  qsave(full_mse, sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_withBez.qs", dep_type, dep_level, box))
}

dep_types <- c("gauss", "logistic", "husler_reiss")
dep_levels <- c("high_wc", "mid_wc", "low_wc", "high", "mid", "low")
boxes <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)
all_combos <- all_combos |> filter(!(dep_types == "gauss" & dep_levels %in% c("low_wc", "mid_wc")),
                                   !(dep_types == "logistic" & dep_levels == "high_wc"),
                                   (!(dep_types == "husler_reiss" & str_detect(dep_levels, "wc"))))

with_progress({
  # Create a progress handler
  p <- progressor(steps = nrow(all_combos))
  
  # Apply the function using apply and update the progress bar
  apply(all_combos, 1, function(row) {
    create_combined_preds(dep_type = row["dep_types"],
                          dep_level = row["dep_levels"],
                          box = row["boxes"])
  })
  p() #update progress bar
})
