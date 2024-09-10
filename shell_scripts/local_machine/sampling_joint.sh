#!/bin/zsh
# model run 
trap '' HUP
basedir="./stan/radial_angular/"
cd ${basedir}
model="bivar_cens_marg_${gauge_name}_angular"

level="mid"
i=1
datafile="../../data/${gauge_name}/${level}_${i}.json"
outbase="csv_fits/${gauge_name}/${level}_${i}_diff_dep"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains for dataset number ${i}, level ${level}, joint likelihood"