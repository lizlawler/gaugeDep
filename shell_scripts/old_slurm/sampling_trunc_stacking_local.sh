#!/bin/zsh
# model run 
trap '' HUP
datafile="../data/${dep_type}/${level}_${i}.json"
basedir="./stan/"
cd ${basedir}
model="bivar_trunc_${gauge_name}"
outbase="csv_fits/stacking/${dep_type}/${gauge_name}/${level}_${i}_trunc_marg"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains"