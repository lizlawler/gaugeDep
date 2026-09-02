#!/bin/zsh
#
# Local (zsh) script: compiles and runs the Stan EGPD-mixture marginal transform model
# to map raw ERC/FWI fire weather data onto exponential margins. This is the Step 0
# prerequisite for the real-data analysis; output goes to samplers/stan/marg_transform/csv_fits/.
# Update the cmdstan stanc path and the data/output file paths for your environment.
#

source ~/.zshrc
conda activate stan
trap '' HUP
stanc_exe="/path/to/your/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="samplers/stan/marg_transform"
# object="samplers/stan/marg_transform/fwi_transform_mix"
object="samplers/stan/marg_transform/erc_fwi_transform_mix_g2"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}

basedir="./samplers/stan/marg_transform/"
cd ${basedir}
model="erc_fwi_transform_mix_g2"

datafile="../../../data/erc_fwi_thomescreek.json"
outbase="csv_fits/transform_thomescreek_g2"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3