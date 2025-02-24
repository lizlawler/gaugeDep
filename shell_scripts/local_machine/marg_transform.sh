#!/bin/zsh
# shell script to test model
#
source  /Users/lizlawler/.zshrc
conda activate stan
trap '' HUP
stanc_exe="/opt/homebrew/Caskroom/miniconda/base/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="samplers/stan/marg_transform"
object="samplers/stan/marg_transform/fwi_transform_mix"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}

basedir="./samplers/stan/marg_transform/"
cd ${basedir}
model="fwi_transform_mix"

datafile="../../../data/erc_fwi.json"
outbase="csv_fits/redstone_fwi_mix_trunc_expo_higher_ub"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3

