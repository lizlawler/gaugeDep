#!/bin/zsh
#
# Local (zsh) script: computes pointwise log-likelihoods for the full sim study by looping
# serially over dep_type, level, gauge, and likelihood. Calls the angular mixture, angular
# star, and radial loglik R scripts. Local serial alternative to the HPC loglik jobs.
#

for dep_type in "gauss" "logistic"; do

	for dep_level in "low" "mid" "high"; do

		export dep_type dep_level

		Rscript --vanilla samplers/nimble/ang_mix_loglik_calc.R ${dep_type} ${dep_level}
		sleep 1
		
		for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do

			export dep_type dep_level gauge

			Rscript --vanilla samplers/rcpp/ang_star_loglik_calc.R ${dep_type} ${dep_level} ${gauge}
			sleep 1

			for likelihood in "trunc" "cens"; do

				export dep_type dep_level gauge likelihood

				Rscript --vanilla samplers/rcpp/radial_loglik_calc.R ${dep_type} ${dep_level} ${gauge} ${likelihood}
				sleep 1

			done

			sleep 1

		done

		sleep 1

	done

	sleep 1

done
