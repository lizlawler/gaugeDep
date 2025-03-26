#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
dep_type="logistic"
level="high"
for batch in {0..4} {10..14} {20..24} {29..34} {39..44}; do
	export dep_type level batch

	sbatch --job-name ${dep_type}_${level}_${batch}_cw_fits_preds \
		--output="./shell_scripts/console_output/%x_%j.txt" \
		shell_scripts/call_cw_logistic_high.sh

	sleep 1
done
