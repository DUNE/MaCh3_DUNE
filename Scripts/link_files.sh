#!/bin/bash
#
# A script to link required files to the proper place i.e. where the sample config files will look for them

MACH3DIR=`pwd`
FILESDIR=/vols/dune/LibanDUNEInputs/
NDGAR_FILESDIR=/vols/dune/jmm224/data

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


if [ ! -d "$MACH3DIR/Inputs/DUNE_NDGAr_files" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_NDGAr_files
fi
ln -sf ${NDGAR_FILESDIR}/NDGAr/rad260/bfield0.5/FastGarSim/FastGArSim_FHC_numu_hA.root Inputs/DUNE_NDGAr_files
if [ ! -d "$MACH3DIR/Inputs/DUNE_NDGAr_files/GenieTrees" ]
then
  mkdir $MACH3DIR/Inputs/DUNE_NDGAr_files/GenieTrees
fi
ln -sf ${NDGAR_FILESDIR}/NDGAr/genie_inputs/*root Inputs/DUNE_NDGAr_files/GenieTrees


ln -sf ${FILESDIR}/DUNE_2023_ND_CAFs_FV_CCINC_Q/*root Inputs/DUNE_ND_CAF_files
