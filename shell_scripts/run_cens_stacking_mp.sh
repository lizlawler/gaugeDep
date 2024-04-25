#!/bin/zsh
#
trap '' HUP
stanc_exe="/data/accounts/lawler/.conda/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/"
# for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"
for gauge_name in "inv_log" "asym_log" "dirichlet" "rectangular"
do
object="stan/bivar_cens_${gauge_name}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for dep_type in "gauss" "logistic"
do
for level in "low" "mid" "high"
do
for i in {1..100}
do
export gauge_name dep_type level i
nohup ./shell_scripts/sampling_cens_stacking.sh > shell_scripts/console_output/stacking/${gauge_name}_${dep_type}_${level}_${i}.txt 2>&1
sleep 3
done
sleep 2
done
sleep 2
done
sleep 2
done