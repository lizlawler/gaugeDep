#!/bin/bash
#
# SLURM job: run the angular mixture (NIMBLE) MCMC sampler over a batch of 5 datasets (sim study).
# Calls: samplers/nimble/ang_mix_mcmc.R
# Env vars in: dep_type, level, batch (start/end dataset index derived from batch)
#

#SBATCH --partition=amilan
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --chdir=/path/to/your/project/gaugeDep/
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@INSTITUTION.EDU

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
export TEMPDIR=${TMPDIR}
mkdir -p $TMPDIR

start_i=$(( batch * 5 + 1 ))
end_i=$(( start_i + 4 ))

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate r_env

Rscript --vanilla samplers/nimble/ang_mix_mcmc.R \
${dep_type} ${level} ${start_i} ${end_i}