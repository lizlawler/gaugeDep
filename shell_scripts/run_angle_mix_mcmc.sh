#!/bin/bash
#
# Submits angular mixture (NIMBLE) MCMC jobs across all sim-study scenarios,
# looping over dep_type, level, and dataset batches. Standalone (no dependency chaining).
#

for dep_type in "gauss" "logistic"; do
  for level in "low" "mid" "high"; do
    for batch in {0..19}; do
      export dep_type level batch
      
      sbatch --job-name ${dep_type}_${level}_${batch}_angle_mix_sampling \
      --output="./shell_scripts/console_output/%x_%j.txt" \
      shell_scripts/call_angle_mix_sampler.sh
      
      sleep 1
    done
      sleep 1
  done
    sleep 1
done
