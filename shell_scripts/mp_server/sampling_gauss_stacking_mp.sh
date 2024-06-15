#!/bin/bash
# model run 

trap '' HUP
basedir="./stan/"
cd ${basedir}
model="bivar_${likelihood}_${threshold}_${gauge_name}"

# run model with 3 chains
for level in "low" "mid" "high" "wc"
do
echo "Starting model runs with for ${level} dependence of AI datasets"
for i in {1..100}
do
datafile="../data/gauss/${level}_${i}.json"
outbase="csv_fits/stacking/gauss/${gauge_name}/${level}_${i}_${likelihood}_${threshold}"
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains for dataset number ${i}"
done
sleep 1
echo "Model has finished running on all 100 datasets for ${level} of AI datasets"
sleep 1
conda activate r_env
Rscript --vanilla extract_params.R \
${dep_type} ${gauge_name} ${likelihood} ${threshold} ${level}
sleep 1
conda activate stan_new
done
sleep 1
echo "${gauge_name}, with ${likelihood} likelihood and ${threshold} threshold has finished running on all datasets"
