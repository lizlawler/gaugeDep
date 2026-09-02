#!/bin/bash
#
# SLURM job: run the BezELS radial MCMC sampler (competitor method, sim study).
# Calls: samplers/bezels/bezels_mcmc.R
# Env vars in: dep_type, level
#

#SBATCH --partition=amilan
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --chdir=/path/to/your/project/gaugeDep/
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@INSTITUTION.EDU

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
export TEMPDIR=${TMPDIR}
mkdir -p $TMPDIR

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate r_env

Rscript --vanilla samplers/bezels/bezels_mcmc.R \
${dep_type} ${level}