#!/bin/bash
#
# A script to link required files to the proper place i.e. where the sample config files will look for them

MACH3DIR=`pwd`
FILESDIR=/vols/dune/ljw20/
FILESDIR1=/vols/dune/nk3717/data/

if [ ! -d "$MACH3DIR/inputs/DUNE_CAF_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_CAF_files
fi
# ln -sf ${FILESDIR}/DUNE_2023_ND_CAFs_FV_CCINC_Q/*root inputs/DUNE_CAF_files
ln -sf ${FILESDIR}/DUNE_2023_FD_CAFs/*root inputs/DUNE_CAF_files


if [ ! -d "$MACH3DIR/inputs/DUNE_spline_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_spline_files
fi

# ln -sf ${FILESDIR}/DUNE_2023_ND_splines/*root inputs/DUNE_spline_files
ln -sf ${FILESDIR}/DUNE_2021_FD_splines/*root inputs/DUNE_spline_files
ln -sf ${FILESDIR}/DUNE_2023_FD_splines/*root inputs/DUNE_spline_files


if [ ! -d "$MACH3DIR/inputs/DUNE_NDGAr_CAF_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_NDGAr_CAF_files
fi
ln -sf ${FILESDIR1}/NDGAr_500kCAFs_2/NDGAr_FHC_ger.root inputs/DUNE_NDGAr_CAF_files


if [ ! -d "$MACH3DIR/inputs/DUNE_NDGAr_spline_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_NDGAr_spline_files
fi
ln -sf ${FILESDIR1}/NDGAr_100kCAFs/SplineOutputs/*root inputs/DUNE_NDGAr_spline_files


if [ ! -d "$MACH3DIR/inputs/DUNE_NDGAr_AnaTrees" ]
then
  mkdir $MACH3DIR/inputs/DUNE_NDGAr_AnaTrees
fi

ln -sf ${FILESDIR}/DUNE_2025_ND_splines/*root inputs/DUNE_ND_spline_files