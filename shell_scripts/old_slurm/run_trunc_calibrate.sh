#!/bin/zsh
#
trap '' HUP
stanc_exe="/opt/homebrew/Caskroom/miniconda/base/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/"
for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet"
do
object="stan/bivar_trunc_${gauge_name}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for level in "low" "mid" "high"
do
for i in {1..100}
do
export gauge_name level i
nohup ./shell_scripts/sampling_trunc_calibrate.sh > shell_scripts/console_output/calibrate/${gauge_name}_${level}_${i}_trunc_ctau.txt 2>&1
sleep 3
done
sleep 2
done
sleep 2
done