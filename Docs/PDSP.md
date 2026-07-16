# PDSP Analysis Workflow

## Running MCMC Fit

```bash
Fit Configs/FitterConfig_PDSP.yaml
```

### Enabling Data Inputs

By default the PDSP fit uses the nominal MC prediction as Asimov data. To fit
data histograms stored in the PDSP input ROOT files, enable:

```yaml
General:
  Data: true
```

or override it at runtime:

```bash
Fit Configs/FitterConfig_PDSP.yaml General:Data:true
```

When `General:Data:true` is set, `SampleHandlerPDSP` looks through each input
ROOT file listed for the sample and loads a histogram ending in `_DataHist`.
The preferred histogram name is the sample title plus `_DataHist`, for example:

```text
PDSP_Abs_DataHist
PDSP_CEx_DataHist
PDSP_Pip_DataHist
```

### Asimov Tunes and MC Scale

Use `XsecAsimovTune` when you want to generate an Asimov data set from the MC
with a known cross-section tune injected before the fit starts. This is useful
for closure tests and generator-consistency checks, where the fake data are
generated with one set of parameter values and the fit then tries to recover
them from the configured prefit values.

The PDSP sample handler also supports a global MC normalization factor:

```yaml
MCGlobalScale: 1.0
```

This scale is read from the sample-handler config and applied as an additional
event weight to every PDSP MC event. Keep it at `1.0` for nominal running, or
change it when the input MC needs a common exposure or normalization
correction.

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

Once you have generated the posterior predictive toy distributions with
PredictivePDSP, you can make fancy plots of them using:

```bash
PredictivePlotting ./Configs/PDSPDiagConfig.yaml PredictiveOutputTest.root
```

### Prior Predictive Distributions

```bash
PredictivePDSP ./Configs/FitterConfig_PDSP.yaml \
  General:OutputFile:PriorPredictiveOutputTest.root \
  Predictive:PriorPredictive:True
```

Finally, we can compare the prior and posterior predictive spectra with the previously used PredictivePlotting macro:

```bash
PredictivePlotting ./Configs/PDSPDiagConfig.yaml PredictiveOutputTest.root PriorPredictiveOutputTest.root
```

## Generator Consistency Checks

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

## Fit Performance Scan

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

To run several scan points concurrently on the local machine:

```bash
./Scripts/run_pdsp_fit_performance_scan.py --jobs 4 --process Abs --process CEx --process Pion
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
