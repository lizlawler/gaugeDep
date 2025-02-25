#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
level="high_wc"
for dep_type in "gauss" "inv_log"; do
  
  mix_mcmc_job_ids=()  # Collect job IDs for angular mixture MCMC runs
  for batch in {0..39}; do
      
    export dep_type level batch
    job_id=$(sbatch --parsable $1 --job-name ${dep_type}_${level}_${batch}_angle_mix_sampling \
                    --output="./shell_scripts/console_output/%x_%j.txt" \
                    shell_scripts/call_angle_mix_sampler.sh)
    mix_mcmc_job_ids+=("$job_id")
    sleep 1
      
  done

  mix_mcmc_dependency_list=$(IFS=,; echo "${mix_mcmc_job_ids[*]}") # Convert job IDs into a comma-separated string, to feed to slurm
  loglik_job_ids=() # Collect job IDs for log likelihood calculation
    
    # Job to calculate pointwise loglikelihood after models have been fit
  job_id=$(sbatch --dependency=afterok:${mix_mcmc_dependency_list} \
                  --parsable $1 --job-name ${dep_type}_${level}_angle_mix_loglik \
                  --output="./shell_scripts/console_output/%x_%j.txt" \
                  shell_scripts/call_angle_mix_loglik_calc.sh)
  loglik_job_ids+=("$job_id") # Collect IDs so model weights can be run automatically once all loglikelihoods have been calculated
  sleep 1 
           
  # Job to extract mean posterior parameters
  sbatch --dependency=afterok:${mix_mcmc_dependency_list} \
         --job-name ${dep_type}_${level}_angle_mix_params \
         --output="./shell_scripts/console_output/%x_%j.txt" \
         shell_scripts/call_ang_mix_params.sh
  sleep 1 
    
  for gauge in "gauss" "logistic" "inv_log" "asym_log" "dirichlet" "rectangular"; do
      
    export dep_type level gauge
    mcmc_job_id=$(sbatch --parsable $1 --job-name ${gauge}_${dep_type}_${level}_angle_star_sampling \
                         --output="./shell_scripts/console_output/%x_%j.txt" \
                         shell_scripts/call_angle_star_sampler.sh)
    sleep 1
      
    # echo "Angular MCMC job ID for ${gauge}, ${dep_type}, ${level}: ${mcmc_job_id}"
    # echo "Submitting loglik job with dependency: afterok:${mcmc_job_id}"
      
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

      # echo "Radial MCMC job ID for ${gauge}, ${dep_type}, ${level}, ${likelihood}: ${mcmc_job_id}"
      # echo "Submitting loglik job with dependency: afterok:${mcmc_job_id}"
        
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
    
  loglik_dependency_list=$(IFS=,; echo "${loglik_job_ids[*]}")
    
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