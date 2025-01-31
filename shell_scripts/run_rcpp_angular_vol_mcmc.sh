#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
for dep_type in "gauss" "logistic"; do
  for level in "low" "mid" "high"; do
    for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do
      export dep_type level gauge_name
      sbatch --job-name ${gauge_name}_${dep_type}_${level}_angle_vol_sampling \
      --output="./shell_scripts/console_output/%x_%j.txt" \
      shell_scripts/call_angle_vol_sampler.sh
      sleep 1
    done
    sleep 1
  done
  sleep 1
done
