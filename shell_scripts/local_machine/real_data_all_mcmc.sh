#!/bin/zsh
#
# Local (zsh) script: runs the full real-data pipeline for both stations (redstone, friendmtn).
# Computes pointwise log-likelihoods for all model components and extracts joint-model BMA
# weights. The MCMC sampling lines are commented out (run separately/beforehand); this script
# focuses on the loglik + weight-extraction stage.
#

for data_type in "redstone" "friendmtn"; do

	export data_type
	# Rscript --vanilla samplers/nimble/ang_mix_mcmc_real_data.R
	Rscript --vanilla samplers/nimble/ang_mix_loglik_calc_real_data.R ${data_type}

	for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do

		export data_type gauge
		# Rscript --vanilla samplers/rcpp/ang_star_mcmc_real_data.R ${gauge}
		# sleep 1

		Rscript --vanilla samplers/rcpp/ang_star_loglik_calc_real_data.R ${data_type} ${gauge}
		sleep 1

		for likelihood in "trunc" "cens"; do

			export data_type gauge likelihood
			# Rscript --vanilla samplers/rcpp/radial_mcmc_real_data.R ${gauge} ${likelihood}
			# sleep 1

			Rscript --vanilla samplers/rcpp/radial_loglik_calc_real_data.R ${data_type} ${gauge} ${likelihood}
			sleep 1

		done

		sleep 1

	done

	for dens in "star" "mix"; do

		for likelihood in "trunc" "cens"; do

			export likelihood dens
			Rscript --vanilla extraction_scripts/extract_weights_joint_real_data.R ${data_type} ${likelihood} ${dens}

		done

		sleep 1

	done

	sleep 1

done
