#!/bin/zsh
#
trap '' HUP

./shell_scripts/run_angle_vol_loglik_calc.sh
./shell_scripts/run_angle_sb_loglik_calc.sh
./shell_scripts/run_radial_loglik_calc.sh