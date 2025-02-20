#!/bin/bash
# shell script to call sbatch
#
# cycle through loop and launch sbatch for every combination
#
dep_type="gauss"
level="high"
likelihood="trunc"
dens="mix"
export dep_type level likelihood dens
sbatch --job-name ${dens}_${dep_type}_${likelihood}_${level}_wts \
       --output="./shell_scripts/console_output/%x_%j.txt" \
       shell_scripts/call_joint_model_wts_extract.sh