#!/bin/bash
#
# Submits BezELS radial MCMC jobs (competitor method) across all sim-study scenarios.
#

for dep_type in "gauss" "logistic" "husler_reiss"; do

	# Set level options based on dep_type
	if [ "$dep_type" == "gauss" ]; then
		levels=("low" "mid" "high" "high_wc")
	elif [ "$dep_type" == "logistic" ]; then
		levels=("low" "mid" "high" "low_wc" "mid_wc")
	elif [ "$dep_type" == "husler_reiss" ]; then
		levels=("low" "mid" "high")
	fi

	for level in "${levels[@]}"; do

		export dep_type level
		sbatch --parsable $1 --job-name ${dep_type}_${level}_radial_bezels_sampling \
			--output="./shell_scripts/console_output/%x_%j.txt" \
			shell_scripts/call_bezels_radial_sampler.sh
		sleep 1

	done

	sleep 1

done
