#!/bin/bash
#
# Submits Campbell-Wadsworth fit-and-predict jobs (competitor method) across all
# sim-study scenarios, in batches of 40 datasets.
#

for dep_type in "gauss" "logistic" "husler_reiss"; do
  for level in "low" "mid" "high"; do
    for batch in {0..4}; do
      export dep_type level batch
      
      sbatch --job-name ${dep_type}_${level}_${batch}_cw_fits_preds \
      --output="./shell_scripts/console_output/%x_%j.txt" \
      shell_scripts/call_cw_fit_and_pred.sh
      
      sleep 1
    done
      sleep 1
  done
    sleep 1
done
