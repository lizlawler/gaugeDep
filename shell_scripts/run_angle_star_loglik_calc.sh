#!/bin/bash
#
# Submits star angular pointwise log-likelihood jobs across all sim-study scenarios and gauges.
#

for dep_type in "gauss" "logistic"; do

  for level in "low" "mid" "high"; do
  
    for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do
    
      export dep_type level gauge
      sbatch --job-name ${gauge}_${dep_type}_${level}_angle_star_loglik \
      --output="./shell_scripts/console_output/%x_%j.txt" \
      shell_scripts/call_angle_star_loglik_calc.sh
      sleep 1
      
    done
    
    sleep 1
    
  done
  
  sleep 1
  
done
