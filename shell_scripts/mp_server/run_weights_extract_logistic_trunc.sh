#!/usr/bin/zsh
# shell script to kick off sampling
#
# cycle through loop and launch sampling for each combination
#
source /data/accounts/lawler/.zshrc

conda activate r_env
trap '' HUP
likelihood="trunc"
dep_type="logistic"
for threshold in "marg" "ctau"
do
for level in "low" "mid"
do
nohup Rscript --vanilla extract_weights.R ${dep_type} ${likelihood} ${threshold} ${level} > shell_scripts/console_output/stacking/${dep_type}_${likelihood}_${threshold}_${level}_weights_mp.txt 2>&1 &
sleep 1
done
sleep 1
done