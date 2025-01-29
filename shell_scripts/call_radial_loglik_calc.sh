#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=csu-general
#SBATCH --chdir=/scratch/alpine/eslawler@colostate.edu/gaugeDep/
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eslawler@colostate.edu

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
export TEMPDIR=${TMPDIR}
mkdir -p $TMPDIR

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate r_env

Rscript --vanilla samplers/rcpp/radial_loglik_calc.R \
${dep_type} ${level} ${gauge_name} ${likelihood}