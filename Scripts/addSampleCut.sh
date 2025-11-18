#!/bin/bash

for entry in $(find /vols/dune/ljw20/DUNE_2023_ND_CAFs/ -name "*.root")
do
  com="root -l -q -b '/vols/t2k/users/ljw20/software/MaCh3_DUNE_ND/MaCh3_DUNE/scripts/addSampleCut.C(\"${entry}\")'"
  eval $com
done
