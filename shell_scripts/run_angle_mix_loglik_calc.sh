#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
for dep_type in "gauss" "logistic"; do

  for level in "low" "mid" "high"; do
  
    export dep_type level
    sbatch --job-name ${dep_type}_${level}_angle_mix_loglik \
    --output="./shell_scripts/console_output/%x_%j.txt" \
    shell_scripts/call_angle_mix_loglik_calc.sh
    sleep 1
    
  done
  
  sleep 1
  
done
