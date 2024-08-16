# MaCh3_DUNE

##################################
# CMAKE #########
##################################

Dependencies

- CMake (version > 3.8). 
- MaCh3 Core tag: DUNECore2024 (To be used until the core version currently being developped gets integrated with MaCh3 DUNE)
- ROOT (currently tested on 6.18)

Building:

~~~~~~~~~~~~~~
$ mkdir MaCh3_DUNE
$ git clone git@github.com:DUNE/MaCh3_DUNE.git MaCh3_DUNE
$ cd MaCh3_DUNE
~~~~~~~~~~~~~~

Now setup some dependencies and then actually build MaCh3_DUNE

~~~~~~~~~~~~~~~
$ source setup.sh # !! Here you need to make sure that ROOTSYS and Cuda libraries are also set !!
$ source setup_dune_env.sh
$ mkdir build;
$ cd build
~~~~~~~~~~~~~~~

Optional flags are described briefly below, options are shown grouped by square brackets and delimited by vertical lines. Default is on the left.

~~~~~~~~~~~~~~
$ cmake .. -DCPU_ONLY=[OFF|ON] -DGPU_ON=[ON|OFF] -DUSE_PROB3=[OFF|ON] -DSINGLE_THREAD_ONLY=[OFF|ON] -DCUDA_SAMPLES=<path_to_cuda>/CentOS/samples
$ make
$ make install
~~~~~~~~~~~~~~

A few notes:
CUDA_SAMPLES not necessary if using CPU_ONLY=ON

If you want to simultaneously develop both the MaCh3 core code and the MaCh3 DUNE code then you can build against a local version of MaCh3 by adding:

~~~~~~~~~~~~~~
$ -DCPM_MaCh3_SOURCE=/path/to/MaCh3/folder
~~~~~~~~~~~~~~

This will overrule the CPMFindPackage command in the CMakeList.txt and will tell CPM to build that instead.

If you have an error about not finding `splineMonolith.h` at the `make install` stage, do: 

~~~~~~~~~~~~~
$ sed -i "s|splineMonolith|SplineMonolith|g" build/_deps/mach3-src/splines/CMakeLists.txt
~~~~~~~~~~~~~ 

###################################
# Event Rates ######
###################################

Once you've got setup you'll then need to setup some symlinks to point to your MC and spline files. You can do this by modifying scripts/link_files.sh script. You'll need to change the FILESDIR variable to point to the relevant folder on your machine. The places these files currently live are listed here:

Imperial College London lx:
~~~~~~~~~~~~~~
/vols/dune/ljw20/
~~~~~~~~~~~~~~

FNAL cluster:
~~~~~~~~~~~~~~
/exp/dune/data/users/lwarsame
~~~~~~~~~~~~~~

ComputeCanada Cedar:
~~~~~~~~~~~~~~
/scratch/liban
~~~~~~~~~~~~~~

NERSC Perlmutter:
~~~~~~~~~~~~~~
/pscratch/sd/l/lwarsame
~~~~~~~~~~~~~~

CVMFS:
~~~~~~~~~~~~~~
/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/MaCh3/inputs/TDR/v2
~~~~~~~~~~~~~~

Current (Feburary 2024) FD event rates using DUNE FD TDR inputs are below (ND is still under-development). These are made using xsec systematics at their prior central value. Oscillation parameter values used here are:

sin2th12 = 0.307

sin2th23 = 0.52

sin2th13 = 0.0218

dm2_32 = 7.53E-5 eV^2

dm2_12 = 2.509E-3 eV^2 

dCP = -1.601 radians

~~~~~~~~~~~~~~~~
Integrals of nominal hists:

FHC_numu unosc:      25941.57467
FHC_numu   osc:      7977.36421
 
FHC_nue unosc:      390.85150
FHC_nue   osc:      1698.28079
 
RHC_numu unosc:      12492.61743
RHC_numu   osc:      4217.78765
 
RHC_nue unosc:      208.31873
RHC_nue   osc:      447.09673
~~~~~~~~~~~~~~~~

