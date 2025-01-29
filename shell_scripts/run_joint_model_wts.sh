#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
for dep_type in "gauss" "logistic"; do
  for level in "low" "mid" "high"; do
    for likelihood in "trunc" "cens"; do
      export dep_type level likelihood
      sbatch --job-name ${likelihood}_${dep_type}_${level}_wts \
      --output="./shell_scripts/console_output/%x_%j.txt" \
      shell_scripts/call_joint_model_wts_extract.sh
      sleep 1
    done
    sleep 1
  done
  sleep 1
done
