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
for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet"
do
for threshold in "marg" "ctau"
do
object="stan/bivar_cens_${threshold}_${gauge_name}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for level in "low" "mid" "high"
do
export gauge_name threshold level
sbatch --job-name ${gauge_name}_${threshold}_${level}_calibration \
--output="./shell_scripts/console_output/calibrate/%x_%j.txt" \
shell_scripts/call_cens_calibrate_sampler.sh
sleep 1
done
sleep 1
done
sleep 1
done
