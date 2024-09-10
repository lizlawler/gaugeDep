#!/bin/bash
# model run 

# change this directory to wherever Stan conda environment lives
source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate stan

basedir="./stan/radial_angular/"
cd ${basedir}
model="bivar_cens_marg_${gauge_name}_mix_betas"

# run model with 3 chains
for i in {1..100}
do
datafile="../../data/${gauge_name}/${level}_${i}.json"
outbase="csv_fits/${gauge_name}/${level}_${i}"
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains for dataset number ${i}"
sleep 1
done
echo "Model has finished running on all 100 datasets"
