#!/bin/bash
#
# A script to link required files to the proper place i.e. where the sample config files will look for them

MACH3DIR=`pwd`
#FILESDIR=/vols/t2k/users/ljw20/data/DUNE_2021/DUNE_2021_splines_tdr_v8
FILESDIR=/vols/dune/ljw20/
FILESDIR1=/vols/dune/nk3717/data/

if [ ! -d "$MACH3DIR/inputs/DUNE_CAF_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_CAF_files
fi
ln -sf ${FILESDIR}/DUNE_2023_FD_CAFs/*root inputs/DUNE_CAF_files


if [ ! -d "$MACH3DIR/inputs/DUNE_spline_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_spline_files
fi

ln -sf ${FILESDIR}/DUNE_2023_FD_splines/*root inputs/DUNE_spline_files


if [ ! -d "$MACH3DIR/inputs/DUNE_NDGAr_CAF_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_NDGAr_CAF_files
fi
ln -sf ${FILESDIR1}/NDGAr_1MCAFs/Outputs/*root inputs/DUNE_NDGAr_CAF_files


if [ ! -d "$MACH3DIR/inputs/DUNE_NDGAr_spline_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_NDGAr_spline_files
fi
ln -sf ${FILESDIR1}/NDGAr_100kCAFs/SplineOutputs/*root inputs/DUNE_NDGAr_spline_files


if [ ! -d "$MACH3DIR/inputs/DUNE_NDGAr_AnaTrees" ]
then
  mkdir $MACH3DIR/inputs/DUNE_NDGAr_AnaTrees
fi
ln -sf ${FILESDIR1}/NDGAr_1MCAFs/AnaTreesOutputs/*root inputs/DUNE_NDGAr_AnaTrees
