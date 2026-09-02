#!/bin/bash
#
# SLURM job: fit Campbell-Wadsworth model + predictions over a batch of 40 datasets (competitor).
# Calls: samplers/campbell_wadsworth/d2fitting.R
# Env vars in: dep_type, level, batch (start/end dataset index derived from batch)
#

#SBATCH --partition=amilan
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --chdir=/path/to/your/project/gaugeDep/
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@INSTITUTION.EDU

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
export TEMPDIR=${TMPDIR}
mkdir -p $TMPDIR

start_i=$(( batch * 40 + 1 ))
end_i=$(( start_i + 39 ))

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate r_env

Rscript --vanilla samplers/campbell_wadsworth/d2fitting.R \
${dep_type} ${level} ${start_i} ${end_i}
