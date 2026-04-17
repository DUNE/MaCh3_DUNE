
source setup_dune_env.sh
export OMP_NUM_THREADS=1

if [ -d build ]; then
    cd build
    cmake .. -DCUDAProb3_ENABLED=OFF -DCUDAProb3Linear_ENABLED=OFF -DBuild_NDGAr=OFF -DDUNE_ANAOBJ_BRANCH="v03_06_00"
    source bin/setup.MaCh3DUNE.sh
    cd ..
else
    echo "build directory not found"
fi