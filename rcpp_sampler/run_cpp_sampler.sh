#!/bin/zsh
#
trap '' HUP
# compile model and link c++ 
# for gauge_name in "gauss" "logistic"
for gauge_name in "logistic"
do
# for lhood in "trunc" "cens"
for lhood in "trunc"
do
for i in {1..100}
do
export gauge_name lhood i
nohup Rscript --vanilla ./rcpp_sampler/sampler_run.R \
${gauge_name} ${lhood} ${i} \
> rcpp_sampler/console_output/${gauge_name}_${lhood}_${i}.txt 2>&1
sleep 3
done
sleep 2
done
sleep 2
done