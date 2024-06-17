#!/usr/bin/zsh
# shell script to kick off sampling
#
# cycle through loop and launch sampling for each combination
#
source /data/accounts/lawler/.zshrc

conda activate r_env
trap '' HUP
char1="cens"
dig1=5
nohup Rscript --vanilla test_rscript.R ${char1} ${dig1} > testing.txt 2>&1 &
