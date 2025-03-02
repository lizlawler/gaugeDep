#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
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

		for likelihood in "trunc" "cens"; do

			for dens in "star" "mix"; do

				export dep_type level likelihood dens
				sbatch --job-name ${dens}_${dep_type}_${likelihood}_${level}_wts \
					--output="./shell_scripts/console_output/%x_%j.txt" \
					shell_scripts/call_joint_model_wts_extract.sh
				sleep 1

			done

			sleep 1

		done

		sleep 1

	done

	sleep 1

done
