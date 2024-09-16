#!/bin/bash
# model run 

# change this directory to wherever Stan conda environment lives
source /curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate stan

basedir="./stan/radial_angular/"
cd ${basedir}
model="bivar_cens_marg_${gauge_name}_mix_betas"

start_i=$(( batch * 20 ))
if [[ $batch -eq 4 ]]; then
    end_i=100  # Ensure the last batch goes through 100
else
    end_i=$(( batch + 19 ))
fi

# run model with 3 chains
for i in $(seq $start_i $end_i) 
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