#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/"
dep_type="logistic"
for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"
do
for likelihood in "trunc" "cens"
do
for threshold in "marg" "ctau"
do
object="stan/bivar_${likelihood}_${threshold}_${gauge_name}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for level in "low" "mid" "high" "wc_mid" "wc_low"
do
export gauge_name likelihood threshold dep_type level
parentjob=$(sbatch --job-name ${gauge_name}_${likelihood}_${threshold}_${dep_type}_${level}_ad_sampling \
--output="./shell_scripts/console_output/stacking/%x_%j.txt" \
shell_scripts/call_${likelihood}_stacking_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${gauge_name}_${likelihood}_${threshold}_${dep_type}_${level}_ad_params \
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
