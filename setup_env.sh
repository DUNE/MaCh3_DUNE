#!/bin/bash
# Full environment setup for MaCh3_DUNE — source this before running any binary.
# Sets up ROOT 6.28 and VDT from CVMFS, then the MaCh3DUNE stack.

ROOT628=/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd
VDT=/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-11.4.1/vdt-0.4.4-kyttjsy6ilqquzjsoi2t2pwiwdoizmnx

export LD_LIBRARY_PATH=${ROOT628}/lib/root:${VDT}/lib:${LD_LIBRARY_PATH}
export ROOTSYS=${ROOT628}

source "$(dirname "$BASH_SOURCE")/build/bin/setup.MaCh3DUNE.sh"
