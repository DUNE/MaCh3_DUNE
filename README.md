# MaCh3_DUNE

## Building MaCh3 DUNE

### Dependencies

- gcc (tested on 12.2.0)
- CMake (tested on 3.27.7)
- ROOT (tested on 6.28.06)

A setup script which pulls cvmfs dependancies is included here:
```bash
source setup_dune_env.sh
```

### Cloning and Building

```bash
mkdir MaCh3_DUNE
git clone git@github.com:DUNE/MaCh3_DUNE.git MaCh3_DUNE
cd MaCh3_DUNE
mkdir build;
cd build
```

Then perform the cmake build command:

```bash
cmake .. -DCUDAProb3_ENABLED=[ON,OFF] -DCUDAProb3Linear_ENABLED=[ON,OFF] -DMaCh3_CORE_BRANCH="v1.4.8"
make install
```

Additional cmake options are available in the MaCh3-Core README

- CUDAProb3 should be used as the default for atmospheric neutrino oscillations
- CUDAProb3Linear should be used as the default for beam oscillations

Then source the installation of MaCh3:
```bash
source build/bin/setup.MaCh3DUNE.sh
```

This sets everything needed, and needs to be re-sourced on each terminal session when using MaCh3 (Along with any dependancies)

## Event Rates

Once you've got setup you'll then need to setup some symlinks to point to your MC and spline files. You can do this by modifying `scripts/link_files.sh` script. You'll need to change the FILESDIR variable to point to the relevant folder on your machine. The places these files currently live are listed here:

Imperial College London lx:
```bash
/vols/dune/ljw20/
```

FNAL cluster:
```bash
/exp/dune/data/users/lwarsame
```

ComputeCanada Cedar:
```bash
/project/rpp-nilic/MaCh3_inputs
```

NERSC Perlmutter:
```bash
/pscratch/sd/l/lwarsame
```

RAL SCARF:
```bash
/work4/ppd/scarf1407
```

CVMFS:
```bash
/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/MaCh3/inputs/TDR/v3
```

Current (Feburary 2024) FD event rates using DUNE FD TDR inputs are below (ND is still under-development). These are made using xsec systematics at their prior central value. Oscillation parameter values used here are:

### Oscillation Parameter Values (NuFIT 4.0 NH)
<div align="center">

|     Parameter     |       Value       |     Unit     |
|:-----------------:|:-----------------:|:------------:|
|     sin²θ₁₂       |       0.310       |      -       |
|     sin²θ₂₃       |       0.582       |      -       |
|     sin²θ₁₃       |       0.0224      |      -       |
|     Δm²₃₂         |    7.39 × 10⁻⁵    |     eV²      |
|     Δm²₁₂         |    2.525 × 10⁻³   |     eV²      |
|     δCP           |      -2.498       |   radians    |

</div>

### Nominal Integrated Event rates

<div align="center">

|       Type        |     Unoscillated    |     Oscillated    |
|:-----------------:|:-------------------:|:-----------------:|
| FHC ν<sub>μ</sub> |     25941.5747      |     8243.9185     |
| FHC ν<sub>e</sub> |      391.5995       |     1756.9128     |
| RHC ν<sub>μ</sub> |     12492.6174      |     4379.4037     |
| RHC ν<sub>e</sub> |      208.8016       |     491.0061      |

</div>
