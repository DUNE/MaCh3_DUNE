#!/bin/bash
#source /home/hep/pjd12/scripts/setup_t2k7.sh
#source setup.sh
folder=$2
nsteps=$3
date=$4
configname=$5
./bin/jointFit2017_5sample_wFGD2_2D $folder/${configname}_${1}.cfg &> $folder/RmuDisapp_2D_${nsteps}steps_${date}_${1}joblog.txt

