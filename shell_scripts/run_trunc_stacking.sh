#!/bin/zsh
#
trap '' HUP
stanc_exe="/opt/homebrew/Caskroom/miniconda/base/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/"
for gauge_name in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"
do
object="stan/bivar_trunc_gamma_${gauge_name}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for dep_type in "independent" "dependent"
do
for level in "low" "mid" "high"
do
for i in {1..100}
do
export gauge_name dep_type level i
nohup ./shell_scripts/sampling_stacking.sh > shell_scripts/console_output/stacking/${gauge_name}_${dep_type}_${level}_${i}.txt 2>&1&
sleep 10
done
sleep 10
done
sleep 10
done
sleep 10
done