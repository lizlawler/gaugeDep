#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
dep_type="husler_reiss"
for level in "mid" " high"; do

	if [ "$level" == "high" ]; then
		batches=({37..40} {78..80} {118..120} {157..160} {198..200})
	elif [ "$level" == "mid" ]; then
		batches=({39..40} {79..80} {119..120} {158..160} {199..200})
	fi

	for batch in "${batches[@]}"; do
		export dep_type level batch

		sbatch --job-name ${dep_type}_${level}_${batch}_cw_fits_preds \
			--output="./shell_scripts/console_output/%x_%j.txt" \
			shell_scripts/call_cw_husler_reiss.sh

		sleep 1
	done

	sleep 1

done
