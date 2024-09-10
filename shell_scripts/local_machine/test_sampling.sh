#!/bin/zsh
# shell script to test model
#
source  /Users/lizlawler/.zshrc
conda activate stan
trap '' HUP
stanc_exe="/opt/homebrew/Caskroom/miniconda/base/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/radial_angular"
gauge_name="gauss"
level="low"
object="stan/radial_angular/bivar_cens_marg_gauss_mix_betas"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
export gauge_name level
basedir="./stan/radial_angular/"
cd ${basedir}
model="bivar_cens_marg_gauss_mix_betas"
datafile="../../data/${gauge_name}/${level}_1.json"
outbase="csv_fits/${gauge_name}/${level}_1_full_model_mix_betas"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv profile_file=${outbase_profile}.csv refresh=25 \
                  num_threads=3



