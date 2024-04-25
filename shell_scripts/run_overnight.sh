#!/bin/zsh
#
trap '' HUP

./shell_scripts/run_trunc_calibrate.sh
./shell_scripts/run_cens_calibrate.sh