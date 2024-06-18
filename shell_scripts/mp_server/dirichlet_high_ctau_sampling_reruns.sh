#!/usr/bin/zsh
# model run 
source /data/accounts/lawler/.zshrc

trap '' HUP
gauge_name="dirichlet"
likelihood="cens"
level="high"
basedir="./stan/"
cd ${basedir}
model="bivar_${likelihood}_${threshold}_${gauge_name}"

# run next on AD datasets
dep_type="logistic"
threshold="ctau"
echo "Starting model runs with for ${level} dependence of AD datasets"
for i in {55..86}
do
datafile="../data/logistic/${level}_${i}.json"
outbase="csv_fits/stacking/logistic/${gauge_name}/${level}_${i}_${likelihood}_${threshold}"
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3

echo "Model has finished running all 3 chains for dataset number ${i}"
done
sleep 1
echo "Model has finished running on all 100 datasets for ${level} of AD datasets"
sleep 1
conda activate r_env
Rscript --vanilla ../extract_params.R \
${dep_type} ${gauge_name} ${likelihood} ${threshold} ${level}
sleep 1
echo "${gauge_name}, with ${likelihood} likelihood and ${threshold} threshold has finished running on all AI and AD datasets"
