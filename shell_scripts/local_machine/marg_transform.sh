#!/bin/zsh
# shell script to test model
#
source  /Users/lizlawler/.zshrc
conda activate stan
trap '' HUP
stanc_exe="/opt/homebrew/Caskroom/miniconda/base/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="samplers/stan/marg_transform"
for family in "g1" "g2" "g3" "g4"
do
object="samplers/stan/marg_transform/fire_transform_${family}"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
og_dir=$(pwd)

basedir="./samplers/stan/marg_transform/"
cd ${basedir}
model="fire_transform_${family}"
# for data in "og" "std"
# do
data="og"
datafile="../../../data/erc_fwi_${data}.json"
outbase="csv_fits/redstone_${data}_${family}"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  init=0.01 \
                  output file=${outbase}.csv \
                  num_threads=3
                  
cd ${og_dir}
done
