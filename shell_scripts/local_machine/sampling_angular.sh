#!/bin/zsh
# model run 
# trap '' HUP
# basedir="./stan/angular/"
# cd ${basedir}
# model="bivar_${gauge_name}_angular"
# 
# for i in {1..100}
# do
# datafile="../../data/angular/${gauge_name}/${level}_${i}.json"
# outbase="csv_fits/${gauge_name}/${level}_${i}"
# 
# # run model with 3 chains
# ./${model} sample num_chains=3 \
#                   data file=${datafile} \
#                   output file=${outbase}.csv \
#                   num_threads=3
#                   
# echo "Model has finished running all 3 chains for dataset number ${i}, level ${level}"
# sleep 1
# done
# echo "Model has finished running on all 100 datasets for ${level} level"

trap '' HUP
basedir="./stan/angular/"
cd ${basedir}
model="bivar_angular_bern"
i=1

for level in "low" "mid" "high" 
do
datafile="../../data/angular/${gauge_name}/${level}_${i}.json"
outbase="csv_fits/${gauge_name}/${level}_${i}_bern"

# run model with 3 chains
./${model} sample num_chains=3 \
                  data file=${datafile} \
                  output file=${outbase}.csv \
                  num_threads=3
                  
echo "Model has finished running all 3 chains for dataset number ${i}, level ${level}, Bernstein density"
sleep 1
done

