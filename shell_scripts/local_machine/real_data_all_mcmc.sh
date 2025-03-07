#!/bin/zsh
# shell script to fit joint model to fire data
#
Rscript --vanilla samplers/nimble/ang_mix_mcmc_real_data.R
Rscript --vanilla samplers/nimble/ang_mix_loglik_calc_real_data.R

for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do

	export gauge
	Rscript --vanilla samplers/rcpp/ang_star_mcmc_real_data.R ${gauge}
	sleep 1

	Rscript --vanilla samplers/rcpp/ang_star_loglik_calc_real_data.R ${gauge}
	sleep 1

	for likelihood in "trunc" "cens"; do

		export gauge likelihood
		Rscript --vanilla samplers/rcpp/radial_mcmc_real_data.R ${gauge} ${likelihood}
		sleep 1

		Rscript --vanilla samplers/rcpp/radial_loglik_calc_real_data.R ${gauge} ${likelihood}
		sleep 1

	done

	sleep 1

done

for dens in "star" "mix"; do

	for likelihood in "trunc" "cens"; do

		export likelihood dens
		Rscript --vanilla extraction_scripts/extract_weights_joint_real_data.R ${likelihood} ${dens}

	done

	sleep 1

done
