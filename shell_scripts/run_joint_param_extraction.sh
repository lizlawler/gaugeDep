#!/bin/bash
#
# Submits posterior parameter extraction jobs (angular mixture, star, radial)
# across all sim-study scenarios. Run after the corresponding MCMC jobs complete.
#

for dep_type in "gauss" "logistic"; do

  for level in "low" "mid" "high"; do

    export dep_type level
    sbatch --job-name ${dep_type}_${level}_ang_mix_params \
    --output="./shell_scripts/console_output/%x_%j.txt" \
    shell_scripts/call_ang_mix_params.sh
    sleep 1
    
    for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do
    
      export dep_type level gauge
      sbatch --job-name ${dep_type}_${level}_${gauge}_ang_star_params \
      --output="./shell_scripts/console_output/%x_%j.txt" \
      shell_scripts/call_ang_star_params.sh
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