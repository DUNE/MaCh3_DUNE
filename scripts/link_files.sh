#!/bin/bash
#
# A script to link required files to the proper place i.e. where the sample config files will look for them
# If running on IC machines this will work out of the box
# If not change FILESDIR to point to where the files live
# Check the README to get information on FILESDIR location according to the cluster used 

MACH3DIR=`pwd`
FILESDIR=/vols/dune/ljw20/

#PATH to /cvmfs where files should live soon!
#FILESDIR=/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/MaCh3/inputs/TDR/v1/

if [ ! -d "$MACH3DIR/inputs/DUNE_CAF_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_CAF_files
fi
ln -sf ${FILESDIR}/DUNE_2023_FD_CAFs/*root inputs/DUNE_CAF_files



if [ ! -d "$MACH3DIR/inputs/DUNE_spline_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_spline_files
fi
ln -sf ${FILESDIR}/DUNE_2021_FD_splines/*root inputs/DUNE_spline_files

if [ ! -d "$MACH3DIR/inputs/DUNE_ND_CAF_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_ND_CAF_files
fi
ln -sf ${FILESDIR}/DUNE_2023_ND_CAFs_FV_CCINC_Q/*root inputs/DUNE_ND_CAF_files



if [ ! -d "$MACH3DIR/inputs/DUNE_ND_spline_files" ]
then
  mkdir $MACH3DIR/inputs/DUNE_ND_spline_files
fi
ln -sf ${FILESDIR}/DUNE_2023_ND_splines/*root inputs/DUNE_ND_spline_files
