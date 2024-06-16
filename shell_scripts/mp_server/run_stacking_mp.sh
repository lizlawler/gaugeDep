#!/usr/bin/zsh
# shell script to kick off sampling
#
# cycle through loop and launch sampling for each combination
#
source /data/accounts/lawler/.zshrc

conda activate stan_new
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan_new/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/"
for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"
do
for likelihood in "trunc" "cens"
do
for threshold in "marg" "ctau"
do
object="stan/bivar_${likelihood}_${threshold}_${gauge_name}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
export gauge_name likelihood threshold
nohup ./shell_scripts/mp_server/sampling_stacking_mp.sh > shell_scripts/console_output/stacking/${gauge_name}_${likelihood}_${threshold}_sampling_mp.txt 2>&1 &
sleep 1
done
sleep 1
done
sleep 1
done
