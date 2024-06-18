#!/usr/bin/zsh
# shell script to kick off sampling
#
# cycle through loop and launch sampling for each combination
#
source /data/accounts/lawler/.zshrc

trap '' HUP
nohup ./shell_scripts/mp_server/asym_log_wclow_ctau_sampling_reruns.sh > shell_scripts/console_output/stacking/asym_log_wclow_ctau_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/asym_log_wclow_marg_sampling_reruns.sh > shell_scripts/console_output/stacking/asym_log_wclow_marg_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/inv_log_wclow_ctau_sampling_reruns.sh > shell_scripts/console_output/stacking/inv_log_wclow_ctau_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/inv_log_wclow_marg_sampling_reruns.sh > shell_scripts/console_output/stacking/inv_log_wclow_marg_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/dirichlet_wclow_ctau_sampling_reruns.sh > shell_scripts/console_output/stacking/dirichlet_wclow_ctau_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/dirichlet_wclow_marg_sampling_reruns.sh > shell_scripts/console_output/stacking/dirichlet_wclow_marg_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/dirichlet_wcmid_ctau_sampling_reruns.sh > shell_scripts/console_output/stacking/dirichlet_wcmid_ctau_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/dirichlet_wcmid_marg_sampling_reruns.sh > shell_scripts/console_output/stacking/dirichlet_wcmid_marg_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/dirichlet_high_ctau_sampling_reruns.sh > shell_scripts/console_output/stacking/dirichlet_high_ctau_sampling_reruns.txt 2>&1 &
sleep 1
nohup ./shell_scripts/mp_server/dirichlet_high_marg_sampling_reruns.sh > shell_scripts/console_output/stacking/dirichlet_high_marg_sampling_reruns.txt 2>&1 &
sleep 1