#!/bin/bash
#
# A script to link required files to the proper place i.e. where the sample config files will look for them

MACH3DIR=`pwd`
FILESDIR=/vols/dune/ljw20/
FILESDIR1=/vols/dune/nk3717/data/

if [ ! -d "$MACH3DIR/Inputs/DUNE_CAF_files" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_CAF_files
fi
ln -sf ${FILESDIR}/DUNE_2023_FD_CAFs/*root Inputs/DUNE_CAF_files


if [ ! -d "$MACH3DIR/Inputs/DUNE_spline_files" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_spline_files
fi

ln -sf ${FILESDIR}/DUNE_2023_FD_splines/*root Inputs/DUNE_spline_files


if [ ! -d "$MACH3DIR/Inputs/DUNE_NDGAr_CAF_files" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_NDGAr_CAF_files
fi
ln -sf ${FILESDIR1}/NDGAr_500kCAFs_2/NDGAr_FHC_ger.root Inputs/DUNE_NDGAr_CAF_files


if [ ! -d "$MACH3DIR/Inputs/DUNE_NDGAr_spline_files" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_NDGAr_spline_files
fi
ln -sf ${FILESDIR1}/NDGAr_100kCAFs/SplineOutputs/*root Inputs/DUNE_NDGAr_spline_files


if [ ! -d "$MACH3DIR/Inputs/DUNE_NDGAr_AnaTrees" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_NDGAr_AnaTrees
fi

if [ ! -d "$MACH3DIR/Inputs/DUNE_ND_spline_files" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_ND_spline_files
fi

ln -sf ${FILESDIR}/DUNE_2025_ND_splines/*root Inputs/DUNE_ND_spline_files

if [ ! -d "$MACH3DIR/Inputs/DUNE_ND_CAF_files" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_ND_CAF_files
fi

ln -sf ${FILESDIR}/DUNE_2023_ND_CAFs_FV_CCINC_Q/*root Inputs/DUNE_ND_CAF_files
