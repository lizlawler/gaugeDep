#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --account=csu-general
#SBATCH --chdir=/scratch/alpine/eslawler@colostate.edu/gaugeDependence/
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=18:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eslawler@colostate.edu

export TMPDIR=/scratch/alpine/$USER/tmp/
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
export TEMPDIR=${TMPDIR}
mkdir -p $TMPDIR

source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate stan

./shell_scripts/sampling_full_model_mixture_betas.sh \
${gauge_name} ${level} ${batch}
