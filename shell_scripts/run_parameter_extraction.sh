#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate lawler

for dep_type in "gauss" "logistic"
do
for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"
do
for likelihood in "trunc" "cens"
do
for threshold in "marg" "ctau"
do
for level in "low" "mid" "high"
do
export dep_type gauge likelihood threshold level
sbatch --job-name ${gauge}_${dep_type}_${likelihood}_${threshold}_${level}_params \
--output="./shell_scripts/console_output/stacking/%x_%j.txt" \
shell_scripts/call_${likelihood}_stacking_params_extraction.sh
sleep 1
done
sleep 1
done
sleep 1
done
sleep 1
done
sleep 1
done
