#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate stan

stanc_exe="/projects/$USER/software/anaconda/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/radial_angular"
for gauge_name in "gauss" "logistic"
do
object="stan/radial_angular/bivar_cens_marg_${gauge_name}_mix_betas"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for level in "low" "mid" "high"
do
export gauge_name level
parentjob=$(sbatch --parsable $1 --job-name ${gauge_name}_${level}_full_lhood_betas_mixture \
--output="./shell_scripts/console_output/%x_%j.txt" \
shell_scripts/call_mixture_betas_sampler.sh)
sleep 1
sbatch --dependency=afterok:${parentjob} \
--job-name ${gauge_name}_${level}_mixture_betas_trace_and_params \
--output="./shell_scripts/console_output/%x_%j.txt" \
shell_scripts/call_mixture_betas_params_extraction.sh
sleep 1
done
sleep 1
done
