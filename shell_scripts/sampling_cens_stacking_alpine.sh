#!/bin/bash
# model run 

# change this directory to wherever Stan conda environment lives
source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate stan

basedir="./stan/"
cd ${basedir}
model="bivar_cens_${threshold}_${gauge_name}"

# run model with 3 chains
for i in {1..100}
do
datafile="../data/${dep_type}/${level}_${i}.json"
outbase="csv_fits/stacking/${dep_type}/${gauge_name}/${level}_${i}_cens_${threshold}"
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains for dataset number ${iter}"
sleep 1
done
echo "Model has finished running on all 100 datasets"
