#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate lawler

dep_type="husler_reiss"
threshold="marg"
for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"
do
for likelihood in "trunc" "cens"
do
for level in "low" "mid" "high"
do
export gauge_name likelihood threshold dep_type level
sbatch --job-name ${gauge_name}_${likelihood}_${threshold}_${dep_type}_${level}_hr_params \
--output="./shell_scripts/console_output/stacking/%x_%j.txt" \
shell_scripts/call_${likelihood}_stacking_params_extraction.sh
sleep 1
done
sleep 1
done
sleep 1
done
