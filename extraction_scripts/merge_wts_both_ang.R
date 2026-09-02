# =============================================================================
# Collects the per-dataset BMA weights produced by extract_weights_joint_both_ang.R
# (one .qs file per dataset) and merges them into a single tidy tibble for
# each combination of dependence type, level, and likelihood. The merged files
# are used for downstream visualisation and summary tables.
#
# Run interactively after all both_ang weight extraction jobs have completed.
# Inputs:    fits_and_weights/wts_joint_model/both_ang/{dep_type}_{dep_level}_{likelihood}_{i}.qs
# Outputs:   fits_and_weights/wts_joint_model/{dep_type}_{likelihood}_{dep_level}_both_ang.qs
# =============================================================================

library(qs)
library(dplyr)
library(tidyr)
library(purrr)

# make_wts_df ----------------------------------------------------------------
# Reads the 200 per-dataset weight files for one dep_type/dep_level/likelihood
# combination, reshapes into a tidy tibble (one row per dataset x model), and
# saves to disk. Short-circuits if the output file already exists.
make_wts_df <- function(dep_type, dep_level, likelihood) {
  wts_file <- sprintf("fits_and_weights/wts_joint_model/%s_%s_%s_both_ang.qs",
                      dep_type, likelihood, dep_level)

  if (file.exists(wts_file)) return(qread(wts_file))

  # 12-element label vector: 6 gauges x 2 angular densities (star, mix)
  gauge_library <- c("gauss", "logistic", "inv_log", "asym_log", "dirichlet", "rectangular")
  ang_lib       <- c("star", "mix")
  model_labels  <- paste(gauge_library, rep(ang_lib, each = 6), sep = "_")

  wts <- map(1:200, function(x) {
    wt_list <- qread(sprintf("fits_and_weights/wts_joint_model/both_ang/%s_%s_%s_%s.qs",
                             dep_type, dep_level, likelihood, x))
    tibble(
      method_ang       = model_labels,
      stacking         = as.numeric(wt_list[[1]]),
      pseudobma_boot   = as.numeric(wt_list[[2]]),
      pseudobma_noboot = as.numeric(wt_list[[3]]),
      dataset          = x
    )
  }) |>
    list_rbind()

  qsave(wts, wts_file)
  print(sprintf("Model weights for %s, %s, %s have been created and saved to disk",
                dep_type, dep_level, likelihood))
}

# Run for all dep_type / dep_level / likelihood combinations
wts_combos <- expand_grid(
  dep   = c("gauss", "logistic"),
  level = c("low", "mid", "high"),
  lhood = c("trunc", "cens")
)

apply(wts_combos, 1, function(row) {
  make_wts_df(dep_type  = row[["dep"]],
              dep_level = row[["level"]],
              likelihood = row[["lhood"]])
})
