#!/bin/zsh
# shell script to kick off sampling
#
# cycle through loop and launch sampling for each combination
#
source  /Users/lizlawler/.zshrc
conda activate stan
trap '' HUP
stanc_exe="/opt/homebrew/Caskroom/miniconda/base/envs/stan/bin/cmdstan/bin/stanc"
# compile model and link c++ 
inc_path="stan/angular"
object="stan/angular/bivar_angular_bern"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
for gauge_name in "gauss" "logistic"
do
export gauge_name
nohup ./shell_scripts/local_machine/sampling_angular.sh > shell_scripts/console_output/angular/${gauge_name}_angular_sampling.txt 2>&1 &
sleep 1
done
