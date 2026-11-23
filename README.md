# MaCh3_DUNE

New here? The [step-by-step guide](#step-by-step-guide-from-setup-to-running-each-app) at the end goes from the environment and the build to running every app (atmospheric samples).

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
cmake .. -DCUDAProb3_ENABLED=[ON,OFF] -DCUDAProb3Linear_ENABLED=[ON,OFF] -DDUNE_ANAOBJ_BRANCH="v03_06_00"
make install
```

Additional cmake options are available in the [MaCh3-Core README](https://github.com/mach3-software/MaCh3?tab=readme-ov-file#other-cmake-options)

- CUDAProb3 should be used as the default for atmospheric neutrino oscillations
- CUDAProb3Linear should be used as the default for beam oscillations

Then source the installation of MaCh3:
```bash
source build/bin/setup.MaCh3DUNE.sh
```

This sets everything needed, and needs to be re-sourced on each terminal session when using MaCh3 (Along with any dependancies)

## Analysis Workflows

- [PDSP analysis workflow](Docs/PDSP.md): fitting, processing MCMC outputs, predictive checks, generator consistency checks, and fit-performance scans.

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
/project/rpp-nilic/MaCh3_Inputs
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

Current (Feburary 2024) FD event rates using DUNE FD TDR Inputs are below (ND is still under-development). These are made using xsec systematics at their prior central value. Oscillation parameter values used here are:

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


## Step-by-step guide: from setup to running each app

This guide walks through the atmospheric analysis (`SampleHandlerAtm`) on the FNAL machines (dunegpvm / EAF). The beam inputs are described in [Event Rates](#event-rates) above.

### 1. Environment (every new shell)

```bash
cd MaCh3_DUNE
source setup_dune_env.sh                 # gcc 12.2.0, ROOT 6.28.12, cmake, vdt from cvmfs
source build/bin/setup.MaCh3DUNE.sh      # after the build: sets MACH3=build/ and the library paths
```

Every path in a config (`Configs/...`, `Inputs/...`) is resolved against `$MACH3`. The build creates the symlinks `build/Configs` and `build/Inputs`, which point at the source tree, so edits to `Configs/` take effect without rebuilding. Apps write their outputs to the current directory.

### 2. Build

```bash
git clone git@github.com:DUNE/MaCh3_DUNE.git MaCh3_DUNE && cd MaCh3_DUNE
git checkout <branch>                    # e.g. dbarrow257/tutorial/ANOAWorkshop
source setup_dune_env.sh
mkdir build && cd build
cmake .. -DCUDAProb3_ENABLED=ON -DCUDAProb3Linear_ENABLED=ON
make -j8 install                         # apps go to build/bin/
cd .. && source build/bin/setup.MaCh3DUNE.sh
```

- CPM downloads the dependencies (MaCh3 core, NuOscillator, CUDAProb3, duneanaobj, yaml-cpp, spdlog, Eigen) into `build/_deps/`.
- After changing code, run `make -j8 install` again in `build/`.
- If the CUDA toolkit isn't found when CMake runs, the GPU code is switched off silently and everything runs on the CPU.

### 3. Inputs

The small atmospheric inputs (binning, Earth model, the 1D detector ratios) are in `Inputs/Atmospherics/`. Link in the large ones. On FNAL they are here:

```bash
SRC=/exp/dune/data/users/barrowd/WorkshopDevelopments/WorkshopTutorial/Inputs/Atmospherics
mkdir -p Inputs/Atmospherics/CAFs Inputs/Atmospherics/Splines
ln -s $SRC/CAFs/AtmSample.eig Inputs/Atmospherics/CAFs/          # MC events, 17 columns (2.1 GB)
ln -s $SRC/DetectorSystematicRatioHistograms_3D.root $SRC/HondaFluxTotalUncertaintySplines.root \
      $SRC/ProdHeightDist_E561_cZ520_N1_hbin2.5.root Inputs/Atmospherics/
ln -s $SRC/Splines/merged_syst.root Inputs/Atmospherics/Splines/ # xsec splines (23.5 GB), only if you switch them on
```

- **MD5 checks:** each input file in `Configs/Samples/SampleHandler_FD_Atm.yaml` has an `...MD5Sum`. A different file stops the app with an MD5 error. For your own file, put `md5sum <file>` in the config.
- **Making the Eigen file yourself:** set `EigenFile: ""` and `FileName` (one CAF) or `FileDirectory` (a folder of CAFs) in the sample config, then run `ConvertToEigen` (step 6).

### 4. Configs

| File | Contents |
|---|---|
| `Configs/EventRates_Atmospherics.yaml` | main config: samples, oscillation parameters, systematics files (`General:Systematics:XsecCovFile`), fixed parameters (`XsecFix`), MCMC settings, plus one section per app (`EventRates`, `PlotSliceHistograms`, `LLHScan`, `Predictive`, `ProcessMCMC`, `Projections`) |
| `Configs/Samples/SampleHandler_FD_Atm.yaml` | atmospheric sample: input files, binning, detector and flux systematics, `FastOptions` |
| `Configs/CovObjs/*.yaml` | parameters (prior, error, bounds, step scale, type): oscillation + Earth model, detector, flux, cross-section |
| `Configs/NuOsc/CUDAProb3_SubSampling_ProdHAvg_EarthSysts.yaml` | oscillation calculation (sub-sampling grid, production-height averaging, Earth model) |

Most apps take the main config, optionally followed by **either** one override YAML **or** command-line overrides:

```bash
./bin/Fit Configs/EventRates_Atmospherics.yaml my_override.yaml           # merged on top; lists replace, they don't append
./bin/Fit Configs/EventRates_Atmospherics.yaml General:OutputFile:Chain.root
```

Command-line overrides of anything under `Systematics` or `Samples` are refused; use an override YAML for those.

**`EventRates`, `PlotSliceHistograms` and `ConvertToEigen` take only the config file** (they stop with "Wrong usage of MaCh3 executable" if given anything more). To change their settings, use a copy of the main config, e.g. `cp Configs/EventRates_Atmospherics.yaml Configs/my_rates.yaml`, and edit it.

By default only the oscillation, detector and flux parameters are loaded, and the cross-section splines are off. **To switch the cross-section splines on:**

1. set `InputSplines: "Inputs/Atmospherics/Splines/merged_syst.root"` in a copy of the sample config;
2. use an override like this one:

```yaml
General:
  OutputFile: Chain_xsec.root
  DUNESamples:
    - Configs/Samples/SampleHandler_FD_Atm_xsec.yaml        # your copy with InputSplines set
  Systematics:
    XsecCovFile:
      - Configs/CovObjs/OscCov_PDG2021_v2_Atmospherics.yaml
      - Configs/CovObjs/DetectorSystematicVariations.yaml
      - Configs/CovObjs/xsec_covariance_Atmospherics.yaml
      - Configs/CovObjs/atmflux_covariance_DUNE.yaml
    XsecFix:
      - ProdH                                                   # no effect on the likelihood, flat prior
      - GENIEReWeight_ICARUS_v2_multisigma_XSecShape_CCMEC         # 0/1 model switches, not 1-sigma dials
      - GENIEReWeight_ICARUS_v2_multisigma_EnergyDependence_CCMEC
```

### 5. Speed options for the atmospheric samples

The fast engine (`FastAtmSplineHandler` and the compact reweight) is on by default and gives the same results as the original code to rounding. The options, with their defaults, are commented in `SampleHandler_FD_Atm.yaml` under `FastOptions`. `SplineEngine: Legacy` together with `FastReweight: false` runs the original code path.

- **Spline cache:** with the splines on, setting `SplineCacheFile` (under `InputFiles`) makes the first run write a binary cache of the knots (3.9 GB instead of 23.5 GB for the 57 dials). Later runs read it instead of the ROOT file, after checking it against `InputSplinesMD5Sum`. You can also make it without building the samples:
  ```bash
  ./bin/MakeAtmSplineCache Inputs/Atmospherics/CAFs/AtmSample.eig Inputs/Atmospherics/Splines/merged_syst.root \
      Configs/CovObjs/xsec_covariance_Atmospherics.yaml Inputs/Atmospherics/Splines/merged_syst.m3spl c35727603ee4f2c58569b906385d31c7 8
  ```
- **Environment overrides** (no config edit): `M3_SPLINE_ENGINE=Fast|Legacy`, `M3_FAST_REWEIGHT=0|1`, `M3_SPLINE_CACHE=<file>|none`, `M3_SPLINE_THREADS=N`, `M3_SPLINE_PRECISION=Double|Float`, `M3_SPLINE_FIXED=Free|Freeze|Drop`.
- **Check that the two paths agree:** `FastCompare` reweights at the nominal point and at reproducible random throws, and writes every bin and −2lnL to a text file:
  ```bash
  M3_SPLINE_ENGINE=Legacy M3_FAST_REWEIGHT=0 FC_OUT=fc_orig.txt ./bin/FastCompare Configs/EventRates_Atmospherics.yaml my_override.yaml
  FC_OUT=fc_fast.txt ./bin/FastCompare Configs/EventRates_Atmospherics.yaml my_override.yaml
  ```
  The other settings are `FC_NTHROWS` (default 20), `FC_SCALE` (throw width in prior σ, 0.5), `FC_REPEAT` (timed reweights per point) and `FC_SEED`.

### 6. Running each app

Run from `build/` after sourcing `build/bin/setup.MaCh3DUNE.sh`. Threads: `export OMP_NUM_THREADS=8`.

| App | What it does | Settings | Output |
|---|---|---|---|
| `EventRates` | builds the samples and fills the event-rate histograms at `General:OscillationParameters`, systematics at prior | `EventRates:PrintIntegral` | `General:OutputFile` (ROOT) and the same name with `.pdf` |
| `PlotSliceHistograms` | draws a multi-dimensional event-rate file as 1D slices | `PlotSliceHistograms` (input file, slice axis, titles) | one PDF per sample (`h<sample>.pdf`) |
| `LikelihoodScan` | 1D −2lnL scan of every free parameter around the Asimov point; 2D scans if `General:2DLLHScan: true` | `LLHScan` (points, ranges, by sample) | ROOT |
| `SigmaVariation` | event-rate histograms with each parameter at −3, −1, 0, +1, +3σ | — | `General:OutputFile` |
| `Fit` | MCMC fit to the Asimov data (or Minuit2 / PSO via `General:FittingAlgorithm`) | `General:MCMC` (`NSteps`, `AutoSave`), `XsecStepScale`, step scales in `Configs/CovObjs` | chain in `General:OutputFile` (tree `posteriors`) |
| `Predictive` | prior predictive (`PriorPredictive: true`) or posterior predictive from a chain (`PosteriorFile`) | `Predictive` | ROOT |
| `Projections` | the samples projected on other variables | `Projections` (`Name`, `VarStrings`, `VarBins`) | ROOT |
| `Variations` | event rates for a list of values of chosen parameters | `Variations` (`Name`, `VarValues`) | ROOT |
| `ConvertToEigen` | atmospheric CAFs → binary Eigen input | sample config `FileName` / `FileDirectory` | `AtmSample.eig` |
| `MakeAtmSplineCache` | spline cache for the fast engine (step 5) | command-line arguments | `.m3spl` |
| `FastCompare` | fast vs original engine, timing and memory (step 5) | `FC_*` environment | text file |

```bash
cd build
./bin/EventRates Configs/EventRates_Atmospherics.yaml                  # writes General:OutputFile = DuneEventRates.root
./bin/PlotSliceHistograms Configs/EventRates_Atmospherics.yaml          # reads PlotSliceHistograms:InputHistogramFileName (DuneEventRates.root)
./bin/LikelihoodScan Configs/EventRates_Atmospherics.yaml General:OutputFile:LLHScan.root
./bin/SigmaVariation Configs/EventRates_Atmospherics.yaml General:OutputFile:SigmaVar.root
./bin/Fit Configs/EventRates_Atmospherics.yaml General:OutputFile:Chain_10k.root   # NSteps from General:MCMC
./bin/Projections Configs/EventRates_Atmospherics.yaml General:OutputFile:Projections.root
./bin/ConvertToEigen Configs/my_caf_config.yaml   # copy of the main config whose DUNESamples sample config has EigenFile: "" and FileName set
```

- **Give each app its own `General:OutputFile`.** The default `DuneEventRates.root` is shared, and `Fit` overwrites (RECREATE) whatever is there.
- **Fit settings:** nested keys work on the command line too, e.g. `./bin/Fit Configs/EventRates_Atmospherics.yaml General:MCMC:NSteps:50000`. Several settings at once go in an override YAML.
- **Variations:** `Configs/Variations_Atmospherics.yaml` is an old sterile-neutrino example; its sample configs no longer exist. To run `Variations`, add a `Variations:` list to your config.
- **Predictive:** it may crash with the atmospheric samples in this core version.

### 7. After the fit: diagnostics and posteriors (MaCh3 core apps, also in `build/bin`)

```bash
./bin/DiagMCMC Chain_10k.root Configs/EventRates_Atmospherics.yaml       # trace, autocorrelation, acceptance -> Chain_10k_MCMC_Diag.root
./bin/PlotMCMCDiag Chain_10k_MCMC_Diag.root diag                          # plots from the diagnostics file
./bin/RHat 10 chain1.root chain2.root                                     # Gelman-Rubin between chains (first argument: thinning)
./bin/CombineMaCh3Chains -o Merged.root chain1.root chain2.root           # merge chains
./bin/ProcessMCMC Configs/EventRates_Atmospherics.yaml Merged.root        # 1D/2D posteriors, credible intervals -> *_drawCorr.root, PDF
./bin/GetPostfitParamPlots Merged_drawCorr.root -o Postfit.pdf            # prior vs posterior per parameter (plotting/ configs, -c to change)
./bin/PlotLLH LLHScan.root                                                # draws the LikelihoodScan output
```

- **Burn-in:** set `ProcessMCMC:BurnInSteps` below the chain length. The default, 80000, is longer than the 10,000-step default fit.
- **Acceptance:** ~0.15–0.35 is healthy. Change `XsecStepScale` or the per-parameter `StepScale` if it is far off.
- **Several chains:** check R̂ < 1.1 before merging.

### 8. Grid jobs (FNAL)

`make_tar.sh` packs the checkout into `../mach3.tar.gz`, and `setup_grid.sh` sources the build's setup scripts on the worker. Send the large inputs as `-f` files (the Eigen file, plus the spline file or the spline cache) rather than inside the tarball.

### 9. Memory and time (full atmospheric sample, 15.25M events, 8 threads)

Measured with the fast engine on the FNAL grid and EAF:

| | memory (peak) | setup | MCMC step |
|---|---|---|---|
| no xsec splines | ~9 GB | — | 0.2–0.4 s |
| 57 xsec splines, fast engine, spline cache | ~13 GB | 35–40 s | 0.3–0.6 s |
| 57 xsec splines, original loader | ~97 GB | ~15 min | 0.9–2.6 s |

The step time is mostly the oscillation calculation. The numbers above include an oscillation speed-up (vectorised production-height averaging in CUDAProb3/NuOscillator) that is not part of this repository. With the stock dependencies, expect that part of each step to be about 1.5× slower.

### 10. Validation (CI)

`Docs/test_EventRates.sh` runs `EventRates` on `CIValidations/CIInputs/CIFitConfig.yaml` and compares the event-rate lines with the stored output (the same check as the pull-request CI).
