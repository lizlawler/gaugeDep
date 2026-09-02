#!/bin/bash
#
# Submits radial (Rcpp) MCMC jobs across all sim-study scenarios, gauges, and both
# likelihoods (trunc, cens). Standalone (no dependency chaining).
#

for dep_type in "gauss" "logistic"; do

  for level in "low" "mid" "high"; do
  
    for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do
    
      for likelihood in "trunc" "cens"; do
      
        export dep_type level gauge likelihood
        sbatch --job-name ${gauge}_${likelihood}_${dep_type}_${level}_radial_sampling \
        --output="./shell_scripts/console_output/%x_%j.txt" \
        shell_scripts/call_${likelihood}_radial_sampler.sh
        sleep 1
        
      done
      
      sleep 1
      
    done
    
    sleep 1
    
  done
  
  sleep 1
  
done
