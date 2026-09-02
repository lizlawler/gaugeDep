# =============================================================================
# Combine level-set (LS) gauge-model predictions with BezELS competitor
# predictions, produce a composite boxplot per scenario, and merge the two
# MSE tables. Run after both imp_samp_preds.R and bezels_preds.R have finished.
#
# Inputs:    figures/is_preds_boxplots/joint/pred_tibbles/{...}.qs
#            figures/is_preds_boxplots/joint/mse_tables/{...}.qs
# Outputs:   figures/is_preds_boxplots/joint/{dep_type}_{dep_level}_{box}_withBez.pdf
#            figures/is_preds_boxplots/joint/plot_objects/{...}_withBez.qs
#            figures/is_preds_boxplots/joint/mse_tables/{...}_withBez.qs
# =============================================================================

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

# Shared ground-truth probability functions
source("sim_study/true_prob_utils.R")

# Load the LS and BezELS prediction tibbles for one scenario, overlay the true
# probability on a log-scale boxplot, and merge the two MSE tables. Box
# dimensions: b1 = symmetric high, b2 = asymmetric, b3 = highly asymmetric.
create_combined_preds <- function(dep_type, dep_level, box) {
  dim1 <- c(10, 12)
  dim2 <- switch(box, b1 = dim1, b2 = c(6, 8), c(2, 4))

  ls_preds     <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s.qs",
                                dep_type, dep_level, box))
  bezels_preds <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_bezels.qs",
                                dep_type, dep_level, box))
  full_preds   <- rbind(ls_preds, bezels_preds)

  true_prob <- get_true_prob(dep_type, dep_level, dim1, dim2)

  # pmax(., eps) keeps zero/near-zero predictions plottable on the log scale.
  plot <- full_preds |>
    ggplot(aes(x = ang_dens, y = pmax(preds, .Machine$double.eps), fill = method)) +
    geom_boxplot() +
    geom_hline(yintercept = pmax(true_prob, .Machine$double.eps),
               col = "darkgrey", linetype = "longdash") +
    scale_y_log10() +
    ggtitle(sprintf("%s, %s, (%s) x (%s)", dep_type, dep_level,
                    paste(dim1, collapse = ","), paste(dim2, collapse = ","))) +
    theme_classic() +
    theme(panel.background  = element_rect(fill = "transparent", color = "transparent"),
          plot.background   = element_rect(fill = "transparent", color = "transparent"),
          axis.text         = element_text(size = rel(1.2)),
          axis.title        = element_text(size = rel(1.2)),
          legend.text       = element_text(size = rel(1.2)),
          legend.title      = element_text(size = rel(1.2)),
          legend.background = element_rect(fill = "transparent", color = "transparent")) +
    scale_fill_grafify(breaks = ~ .x[!is.na(.x)], palette = "r4", ColSeq = FALSE) +
    labs(fill = "BMA method") +
    scale_x_discrete(labels = c("BezELS", "Cens., Mix", "Cens., Star", "Trunc., Mix", "Trunc., Star")) +
    xlab("Likelihood, Angular Density") + ylab("Prediction probabilities")

  ggsave(sprintf("figures/is_preds_boxplots/joint/%s_%s_%s_withBez.pdf", dep_type, dep_level, box),
         plot = plot, bg = "transparent", width = 8, height = 7, dpi = 320)
  qsave(plot, sprintf("figures/is_preds_boxplots/joint/plot_objects/%s_%s_%s_withBez.qs",
                      dep_type, dep_level, box))
  print(sprintf("Joint boxplot for %s %s, in box %s has been saved", dep_level, dep_type, box))

  # Merge and re-sort the MSE tables. (Sort key is rmse_norm, the column the
  # MSE tables actually contain.)
  ls_mse   <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s.qs",
                            dep_type, dep_level, box))
  bez_mse  <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_bezels.qs",
                            dep_type, dep_level, box))
  full_mse <- rbind(ls_mse, bez_mse) |> arrange(rmse_norm)
  qsave(full_mse, sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_withBez.qs",
                          dep_type, dep_level, box))
}

# Scenario grid. The _wc ("well-classified") levels exist only for certain
# dep_types, so the invalid combinations are filtered out.
dep_types  <- c("gauss", "logistic", "husler_reiss")
dep_levels <- c("high_wc", "mid_wc", "low_wc", "high", "mid", "low")
boxes      <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes) |>
  filter(!(dep_types == "gauss"        & dep_levels %in% c("low_wc", "mid_wc")),
         !(dep_types == "logistic"     & dep_levels == "high_wc"),
         !(dep_types == "husler_reiss" & str_detect(dep_levels, "wc")))

with_progress({
  p <- progressor(steps = nrow(all_combos))
  apply(all_combos, 1, function(row) {
    create_combined_preds(dep_type  = row["dep_types"],
                          dep_level = row["dep_levels"],
                          box       = row["boxes"])
    p()
  })
})
