#!/bin/zsh
# model run 
trap '' HUP
datafile="../data/${gauge_name}/${level}_${i}.json"
basedir="./stan/"
cd ${basedir}
model="bivar_trunc_gamma_${gauge_name}"
outbase="csv_fits/calibrate/${gauge_name}/${level}_${i}"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains"