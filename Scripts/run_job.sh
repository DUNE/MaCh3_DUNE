#!/bin/bash
#source /home/hep/pjd12/scripts/setup_t2k7.sh
#source setup.sh
folder=$2
nsteps=$3
date=$4
./bin/jointFit2015_wFGD2 $folder/kirstyconfig_1D_${1}.cfg &> $folder/RmuDisapp_1D_${nsteps}steps_${date}_${1}joblog.txt

