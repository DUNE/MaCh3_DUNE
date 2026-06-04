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
source setup.sh 
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
