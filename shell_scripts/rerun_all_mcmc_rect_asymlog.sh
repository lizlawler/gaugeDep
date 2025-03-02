#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
for dep_type in "gauss" "logistic"; do

	# Set level options based on dep_type
	if [ "$dep_type" == "gauss" ]; then
		levels=("low" "mid" "high" "high_wc")
	elif [ "$dep_type" == "logistic" ]; then
		levels=("low" "mid" "high" "low_wc" "mid_wc")
	fi

	for level in "${levels[@]}"; do

		loglik_job_ids=() # Collect job IDs for log likelihood calculation

		for gauge in "asym_log" "rectangular"; do

			export dep_type level gauge
			mcmc_job_id=$(sbatch --parsable $1 --job-name ${gauge}_${dep_type}_${level}_angle_star_sampling \
				--output="./shell_scripts/console_output/%x_%j.txt" \
				shell_scripts/call_angle_star_sampler.sh)
			sleep 1

			job_id=$(sbatch --dependency=afterok:${mcmc_job_id} \
				--parsable $1 --job-name ${gauge}_${dep_type}_${level}_angle_star_loglik \
				--output="./shell_scripts/console_output/%x_%j.txt" \
				shell_scripts/call_angle_star_loglik_calc.sh)
			loglik_job_ids+=("$job_id")
			sleep 1

			sbatch --dependency=afterok:${mcmc_job_id} \
				--job-name ${dep_type}_${level}_${gauge}_ang_star_params \
				--output="./shell_scripts/console_output/%x_%j.txt" \
				shell_scripts/call_ang_star_params.sh
			sleep 1

			for likelihood in "trunc" "cens"; do

				export dep_type level gauge likelihood
				mcmc_job_id=$(sbatch --parsable $1 --job-name ${gauge}_${likelihood}_${dep_type}_${level}_radial_sampling \
					--output="./shell_scripts/console_output/%x_%j.txt" \
					shell_scripts/call_${likelihood}_radial_sampler.sh)
				sleep 1

				job_id=$(sbatch --dependency=afterok:${mcmc_job_id} \
					--parsable $1 --job-name ${gauge}_${likelihood}_${dep_type}_${level}_radial_loglik \
					--output="./shell_scripts/console_output/%x_%j.txt" \
					shell_scripts/call_radial_loglik_calc.sh)
				loglik_job_ids+=("$job_id")
				sleep 1

				sbatch --dependency=afterok:${mcmc_job_id} \
					--job-name ${gauge}_${likelihood}_${dep_type}_${level}_radial_params \
					--output="./shell_scripts/console_output/%x_%j.txt" \
					shell_scripts/call_radial_params.sh
				sleep 1

			done

			sleep 1

		done

		loglik_dependency_list=$(
			IFS=,
			echo "${loglik_job_ids[*]}"
		)

		for dens in "star" "mix"; do

			for likelihood in "trunc" "cens"; do

				export dep_type level likelihoood dens
				sbatch --dependency=afterok:${loglik_dependency_list} \
					--job-name ${dens}_${dep_type}_${level}_wts \
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
