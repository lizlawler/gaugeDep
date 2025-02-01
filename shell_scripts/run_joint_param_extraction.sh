#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
for dep_type in "gauss" "logistic"; do
  for level in "low" "mid" "high"; do
    export dep_type level
    sbatch --job-name ${dep_type}_${level}_sb_params \
    --output="./shell_scripts/console_output/%x_%j.txt" \
    shell_scripts/call_sb_params.sh
    sleep 1
    for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do
      export dep_type level gauge
      sbatch --job-name ${dep_type}_${level}_${gauge}_vol_params \
      --output="./shell_scripts/console_output/%x_%j.txt" \
      shell_scripts/call_vol_params.sh
      sleep 1
      for likelihood in "trunc" "cens"; do
        export dep_type level gauge likelihood
        sbatch --job-name ${gauge}_${likelihood}_${dep_type}_${level}_radial_params \
        --output="./shell_scripts/console_output/%x_%j.txt" \
        shell_scripts/call_radial_params.sh
        sleep 1
      done
    sleep 1
    done
    sleep 1
  done
  sleep 1
done
  
