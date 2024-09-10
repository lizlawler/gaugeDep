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
# for gauge_name in "gauss" "logistic"
# do
gauge_name="logistic"
inc_path="stan/radial_angular"
object="stan/radial_angular/bivar_cens_marg_${gauge_name}_angular"
${stanc_exe} ${object}.stan --include-paths=${inc_path}
cmdstan_model ${object}
export gauge_name
nohup ./shell_scripts/local_machine/sampling_joint.sh > shell_scripts/console_output/joint/${gauge_name}_joint_sampling.txt 2>&1 &
# sleep 1
# done
