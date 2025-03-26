#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
dep_type="gauss"
# level="high"
level="mid"
for batch in {36..40} {77..80}; do
	export dep_type level batch

	sbatch --job-name ${dep_type}_${level}_${batch}_cw_fits_preds \
		--output="./shell_scripts/console_output/%x_%j.txt" \
		shell_scripts/call_cw_gauss.sh

	sleep 1
done