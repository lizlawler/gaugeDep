#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=csu-general
#SBATCH --chdir=/scratch/alpine/eslawler@colostate.edu/gaugeDep/
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eslawler@colostate.edu

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
export TEMPDIR=${TMPDIR}
mkdir -p $TMPDIR

start_i=$(( batch * 4 + 1 ))
end_i=$(( start_i + 3 ))

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate r_env

Rscript --vanilla samplers/nimble/angular_sb_mcmc.R \
${dep_type} ${level} ${start_i} ${end_i}