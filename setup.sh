
source setup_dune_env.sh

ROOT628=/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd
VDT=/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-11.4.1/vdt-0.4.4-kyttjsy6ilqquzjsoi2t2pwiwdoizmnx

export LD_LIBRARY_PATH=${ROOT628}/lib/root:${VDT}/lib:${LD_LIBRARY_PATH}
export ROOTSYS=${ROOT628}
export OMP_NUM_THREADS=1

if [ -d build ]; then
    cd build
    cmake .. -DCUDAProb3_ENABLED=OFF -DCUDAProb3Linear_ENABLED=OFF -DBuild_NDGAr=OFF -DDUNE_ANAOBJ_BRANCH="v03_06_00"
    source bin/setup.MaCh3DUNE.sh
    cd ..
else
    echo "build directory not found"
fi