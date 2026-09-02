#!/bin/bash
#
# SLURM job: compute joint-model BMA weights (stacking / pseudo-BMA) via loo (sim study).
# Calls: extraction_scripts/extract_weights_joint.R
# Env vars in: dep_type, level, likelihood, dens
#

#SBATCH --partition=amilan
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --chdir=/path/to/your/project/gaugeDep/
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@INSTITUTION.EDU

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
export TEMPDIR=${TMPDIR}
mkdir -p $TMPDIR

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate r_env

Rscript --vanilla extraction_scripts/extract_weights_joint.R \
${dep_type} ${level} ${likelihood} ${dens}