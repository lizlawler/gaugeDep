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
for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"
do
for threshold in "marg" "ctau"
do
object="stan/bivar_trunc_${threshold}_${gauge_name}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for dep_type in "gauss" "logistic"
do
for level in "low" "mid" "high"
do
export gauge_name threshold dep_type level
sbatch --job-name ${gauge_name}_${threshold}_${dep_type}_${level}_trunc \
--output="./shell_scripts/console_output/stacking/%x_%j.txt" \
shell_scripts/call_trunc_stacking_sampler.sh
sleep 1
done
sleep 1
done
sleep 1
done
sleep 1
done
