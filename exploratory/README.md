# Exploratory scripts

These scripts are kept for local reference only and are **excluded from the
public repository** via `.gitignore`. They are exploratory or superseded work
from earlier stages (prelim/defense presentations and analysis iterations) and
are **not** part of the reproducible manuscript pipeline. Some may not run
as-is against the current code/data.

- `bivar_data_viz.R` — early bivariate data-illustration figures (superseded by
  `viz_scripts/viz_copulas.R`-style plots; produced prelim-presentation figures).
- `viz_copulas.R` — exploratory copula/gauge visualizations (not in the paper).
- `viz_pred_task.R` — exploratory prediction-task visualizations (useful but not
  the final manuscript figure).
- `viz_exceed_quant_reg.R` — near-duplicate of
  `viz_scripts/viz_conditional_exceedances_true_v_gamma.R`; only the latter
  produces the manuscript figure.
- `stacking_wts_boxplots.R` — exploratory BMA stacking-weight boxplots (not in
  the paper).
- `imp_samp_preds_real_data.R` — superseded by
  `data_analyses/imp_samp_real_data_all_iter.R`. The newer script runs
  predictions over all posterior iterations and also revises the method:
  station-specific prediction boxes read from `pred_boxes.qs`, an importance
  proposal scaled to each box, and a truncation correction for the positive
  quadrant. Results therefore differ from this script even at the median.
