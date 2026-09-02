# =============================================================================
# Combine predictions from every model variant -- LS with each fixed angular
# density, LS with both angular densities jointly ("Both"), BezELS, and C-W --
# into one composite boxplot per scenario, and merge their MSE tables. Run
# after imp_samp_preds.R, imp_samp_preds_both_ang.R, bezels_preds.R, and
# reshape_cw_preds.R have all completed.
#
# Inputs:    figures/is_preds_boxplots/joint/pred_tibbles/{...}.qs
#            figures/is_preds_boxplots/joint/mse_tables/{...}.qs
# Outputs:   figures/is_preds_boxplots/joint/{dep_type}_{dep_level}_{box}_{all|all_winsor}_both_ang.pdf
#            figures/is_preds_boxplots/joint/plot_objects/{...}_both_ang.qs
#            figures/is_preds_boxplots/joint/mse_tables/{...}_all_both_ang.qs
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

# Load the LS, both-angular, BezELS, and C-W prediction tibbles for one
# scenario; optionally winsorise the predictions at the 2.5/97.5 percentiles;
# draw a composite boxplot with the true probability overlaid; and merge the
# MSE tables. Box dimensions: b1 = symmetric high, b2 = asymmetric,
# b3 = highly asymmetric.
create_combined_preds <- function(dep_type, dep_level, box, plot = TRUE, winsor = FALSE) {
  dim1 <- c(10, 12)
  dim2 <- switch(box, b1 = dim1, b2 = c(6, 8), c(2, 4))

  ls_preds       <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s.qs",
                                  dep_type, dep_level, box))
  both_ang_preds <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/both_ang/%s_%s_cens_%s.qs",
                                  dep_type, dep_level, box)) |>
    mutate(ang_dens = "Cens., Both", .before = "method") |>
    select(-model)
  bezels_preds   <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_bezels.qs",
                                  dep_type, dep_level, box))
  cw_preds       <- qread(sprintf("figures/is_preds_boxplots/joint/pred_tibbles/%s_%s_%s_cw.qs",
                                  dep_type, dep_level, box))

  # Standardise the ang_dens labels and fix their display order.
  full_preds <- rbind(ls_preds, both_ang_preds, bezels_preds, cw_preds) |>
    mutate(ang_dens = case_when(
      ang_dens == "cens, mix"   ~ "Cens., Mix",
      ang_dens == "cens, star"  ~ "Cens., Star",
      ang_dens == "trunc, mix"  ~ "Trunc., Mix",
      ang_dens == "trunc, star" ~ "Trunc., Star",
      ang_dens == "BezELS"      ~ "(BezELS)",
      ang_dens == "C-W"         ~ "(C-W)",
      .default = ang_dens
    )) |>
    mutate(ang_dens = factor(ang_dens,
                             levels = c("Cens., Both",
                                        "Cens., Mix", "Cens., Star",
                                        "Trunc., Mix", "Trunc., Star",
                                        "(BezELS)", "(C-W)")))

  true_prob <- get_true_prob(dep_type, dep_level, dim1, dim2)

  if (plot) {
    fileend <- if (winsor) {
      # Winsorise per group to keep extreme outliers from dominating the plot.
      winsor_limits <- full_preds |>
        group_by(ang_dens, method) |>
        summarize(lower = quantile(preds, 0.025), upper = quantile(preds, 0.975),
                  .groups = "drop")
      full_preds <- full_preds |>
        left_join(winsor_limits, by = c("ang_dens", "method")) |>
        mutate(preds = case_when(preds <= lower ~ lower,
                                 preds >= upper ~ upper,
                                 .default = preds))
      "all_winsor"
    } else { "all" }

    # Numeric x positions so the BMA boxes (1-5) and the singleton BezELS/C-W
    # boxes (offset at 5.75/6.25) share consistent spacing.
    x_map <- tibble(
      ang_dens = c("Cens., Both", "Cens., Mix", "Cens., Star",
                   "Trunc., Mix", "Trunc., Star", "(BezELS)", "(C-W)"),
      x_pos    = c(1, 2, 3, 4, 5, 5.75, 6.25)
    )
    ls_data    <- full_preds |> filter(!is.na(method)) |> left_join(x_map, by = "ang_dens")
    other_data <- full_preds |> filter( is.na(method)) |> left_join(x_map, by = "ang_dens")

    plot_obj <- ggplot() +
      geom_boxplot(data = ls_data,
                   aes(x = x_pos, y = pmax(preds, .Machine$double.eps),
                       fill = method, group = interaction(x_pos, method)),
                   width = 0.75) +
      geom_boxplot(data = other_data,
                   aes(x = x_pos, y = pmax(preds, .Machine$double.eps), group = x_pos),
                   fill = "transparent", color = "black", width = 0.25) +
      geom_hline(yintercept = pmax(true_prob, .Machine$double.eps),
                 col = "darkgrey", linetype = "longdash") +
      scale_x_continuous(breaks = x_map$x_pos, labels = x_map$ang_dens) +
      scale_y_log10() +
      scale_fill_grafify(breaks = ~ .x[!is.na(.x)], palette = "r4", ColSeq = FALSE) +
      labs(fill = "BMA method",
           x    = "Likelihood, Angular Density (Model)",
           y    = "Prediction probabilities") +
      theme_classic() +
      theme(panel.background  = element_rect(fill = "transparent", color = "transparent"),
            plot.background   = element_rect(fill = "transparent", color = "transparent"),
            axis.text         = element_text(size = rel(1.2)),
            axis.title        = element_text(size = rel(1.2)),
            legend.text       = element_text(size = rel(1.2)),
            legend.title      = element_text(size = rel(1.2)),
            legend.background = element_rect(fill = "transparent", color = "transparent"))

    ggsave(sprintf("figures/is_preds_boxplots/joint/%s_%s_%s_%s_both_ang.pdf",
                   dep_type, dep_level, box, fileend),
           plot = plot_obj, bg = "transparent", width = 12, height = 7, dpi = 320)
    qsave(plot_obj, sprintf("figures/is_preds_boxplots/joint/plot_objects/%s_%s_%s_%s_both_ang.qs",
                            dep_type, dep_level, box, fileend))
    print(sprintf("Joint boxplot for %s %s, in box %s has been saved", dep_level, dep_type, box))
  }

  # Merge and re-sort the MSE tables (sort key rmse_norm, as the tables contain).
  ls_mse       <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s.qs",
                                dep_type, dep_level, box))
  both_ang_mse <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/both_ang/%s_%s_cens_%s.qs",
                                dep_type, dep_level, box)) |>
    mutate(ang_dens = "cens, Both", .before = "mse")
  bez_mse      <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_bezels.qs",
                                dep_type, dep_level, box))
  cw_mse       <- qread(sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_cw.qs",
                                dep_type, dep_level, box))
  full_mse <- rbind(ls_mse, both_ang_mse, bez_mse, cw_mse) |> arrange(rmse_norm)
  qsave(full_mse, sprintf("figures/is_preds_boxplots/joint/mse_tables/%s_%s_%s_all_both_ang.qs",
                          dep_type, dep_level, box))
}

# The both-angular predictions exist only for gauss and logistic.
dep_types  <- c("gauss", "logistic")
dep_levels <- c("high", "mid", "low")
boxes      <- c("b1", "b2", "b3")
all_combos <- expand_grid(dep_types, dep_levels, boxes)

with_progress({
  p <- progressor(steps = nrow(all_combos))
  apply(all_combos, 1, function(row) {
    create_combined_preds(dep_type  = row["dep_types"],
                          dep_level = row["dep_levels"],
                          box       = row["boxes"],
                          plot      = TRUE,
                          winsor    = TRUE)
    p()
  })
})
