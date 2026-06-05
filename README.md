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

## Running MCMC fit
```bash
Fit Configs/FitterConfig_PDSP.yaml
```

## Processing MCMC Outputs
```bash
ProcessMCMC ./Configs/PDSPDiagConfig.yaml Test.root
```

## Posterior Predictive Analysis
Once you run MCMC you can produce these toy distributions using following command:
```bash
PredictivePDSP Configs/FitterConfig_PDSP.yaml General:OutputFile:PredictiveOutputTest.root
```

### Plotting Posterior Predictive Distributions
Once you have generated the posterior predictive toy distributions with PredictivePDSP, you can make fancy plots of them using:
```bash
PredictivePlotting ./Configs/PDSPDiagConfig.yaml PredictiveOutputTest.root
```

### Prior Predictive Distributions
```bash
PredictivePDSP ./Configs/FitterConfig_PDSP.yaml General:OutputFile:PriorPredictiveOutputTest.root Predictive:PriorPredictive:True
```

Finally, we can compare the prior and posterior predictive spectra with the previously used PredictivePlotting macro:
```bash
PredictivePlotting ./Configs/PDSPDiagConfig.yaml PredictiveOutputTest.root PriorPredictiveOutputTest.root
```

## PDSP Generator consistency checks
To run process-level varied checks, use:
```bash
./Scripts/run_pdsp_generator_consistency.py --set Abs=1.2 --set CEx=0.8
```
By default this runs the full result chain for each setting: `Fit`,
`ProcessMCMC`, posterior `PredictivePDSP`, prior `PredictivePDSP`, and
`PredictivePlotting`. It writes copied configs, ROOT outputs, logs, per-case
manifests, and a summary table under `PDSPGeneratorConsistency/`. The source
`Configs/CovObjs/PDSPFitModel.yaml` is not modified.
The overlay plots from `PredictivePlotting` are written in each case directory,
for example `PDSPGeneratorConsistency/Abs_generator_1.2/Overlay_Predictive.pdf`.

To also run a nominal fake-data check with `Generator = 1` for every PDSP
systematic:
```bash
./Scripts/run_pdsp_generator_consistency.py --nominal --set Abs=1.2
```

To only generate the fit chain and skip the predictive/plotting steps:
```bash
./Scripts/run_pdsp_generator_consistency.py --workflow fit --set Abs=1.2
```
List available process and parameter names with:
```bash
./Scripts/run_pdsp_generator_consistency.py --list
```
You can also vary one systematic exactly:
```bash
./Scripts/run_pdsp_generator_consistency.py --parameter Abs_TrueEBin_0=1.5
```

## PDSP fit performance scan
To scan injected cross-section normalisations and estimate where fit recovery
starts to degrade, use:
```bash
./Scripts/run_pdsp_fit_performance_scan.py --process Abs --process CEx --process Pion
```
By default this scans `Generator = 1.1, 1.2, ..., 2.0` for each requested
process and runs the fit stage only. It then compares each varied parameter's
posterior mean against the injected value after burn-in and writes:
```text
PDSPFitPerformanceScan/performance_summary.csv
PDSPFitPerformanceScan/safe_region_summary.csv
```
The default degradation criterion is `abs((posterior_mean - injected) /
injected) > 0.10`. Change this with `--tolerance`, for example:
```bash
./Scripts/run_pdsp_fit_performance_scan.py --process Abs --tolerance 0.05
```
To include the nominal `Generator = 1.0` point explicitly:
```bash
./Scripts/run_pdsp_fit_performance_scan.py --include-nominal --process Abs
```
To run the full predictive plotting chain at every scan point:
```bash
./Scripts/run_pdsp_fit_performance_scan.py --workflow full --process Abs
```
To visualise degradation versus injected normalisation after the scan:
```bash
./Scripts/plot_pdsp_fit_performance_scan.py
```
This reads `PDSPFitPerformanceScan/performance_summary.csv` and writes one plot
per scanned process or parameter under `PDSPFitPerformanceScan/plots/`, with
injected `Generator` normalisation on the x-axis and relative discrepancy on the
y-axis. Vertical error bars show `posterior_rms / injected`, i.e. the posterior
width on the same relative scale. If the top-level summary CSV is missing or
does not contain the requested target, the plotter rebuilds it from the
per-point `*_recovery.csv` files.

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
