# BeamOffAxis Analysis

## Inputs

### Sample yml

SampleHandler: Contains multiple *Samples* -- we will use these samples as detector positions
-- Each *Sample* contains multiple 'subsamples' (which correspond to oscillation channels, we will just use a single subsample per sample.)

Can specify POT per *Sample*

#### Binning

Uniform: True for now

-- short description of VarBins options

Global bin number across *Samples* -- details don't matter to analyser except for covariance matrix.

### XSec Splines

Generate splines by passing spline yml

```bash
SplineMaker Configs/BeamOffAxis/EventRate.yml Configs/BeamOffAxis/XSecSysts.yml
```

Don't need to comment out Systematics.XsecCov in EventRate.yml, it is ignored by SplineMaker

This outputs one spline file per *Sample* defined in the SampleHandler into current working dir.

  e.g.: `OnAxis_numuCC_numode_splines.root`, `OffAxis8m_numuCC_numode_splines.root`

Commented out a bunch of parameters that had missing branches/no response.

### ND Detector Covariance Matrix

Make Det cov like:

```bash
MakeCovarianceMatrix Configs/BeamOffAxis/EventRate.yml Configs/BeamOffAxis/NDDetSysts.yml 1000 bla.root
```

will do 1000 throws and write to bla.root. Covariance matrix is called: `NDCovMatrix`.
Assumes and only works with a single configured SampleHandler (but multiple *Samples*).

### Template Parameter definition

New python script for generating:

```bash
python3 ../Samples/BeamOffAxis/scripts/gen_template_param_yml.py > ../Configs/BeamOffAxis/TemplateParams.yml
```

## Anaysis Details

Now Samples/SampleHandlerBeamOffAxis.{h,cpp} has minimal implementation, much has been moved
to smaller files in Samples/BeamOffAxis/

### Event Reading

See Samples/BeamOffAxis/EventInfo.h for event structure
See Samples/BeamOffAxis/ReadEvents.h for code that reads events from CAF and calculates some
composite properties.

Other analysis-specific properties are calculated in Sample/SampleHandlerBeamOffAxis.cpp:SetupExperimentMC

### Projections

Samples/BeamOffAxis/Projections.h for ReturnKinematicParameter

### Variational Systematics

Based on multi-parameter functional variations

-- FinalizeShifts calculates ERec in Samples/BeamOffAxis/Systematics.cpp:CalculateVariedCompositeQuantities after all the individual shift parameters have run.

### Likelihood

Likelihood defined in Samples/SampleHandlerBeamOffAxis.h:GetLikelihood

-- Luke needs to work out how to do regularisation.
