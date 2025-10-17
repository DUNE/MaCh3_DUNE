# Flux Systematic Handling

## Input Provenance

The flux ratios used here were taken from CAFAna/prism in October 2025. The 
inputs we use here were synthesised by combining two sets of histograms, an 
older set ([`flux_shifts_OffAxis.root`](https://github.com/DUNE/lblpwgtools/blob/b9fd2dc0acca52a6c20a246cb3ceea313c0a4529/CAFAna/Systs/flux_shifts_OffAxis.root)) 
containing the flux ratios from throws of the hadron production model (via 
the package PPFX) and an updated set 
([`flux_shifts_OffAxis2023.root`](https://github.com/DUNE/lblpwgtools/blob/b9fd2dc0acca52a6c20a246cb3ceea313c0a4529/CAFAna/Systs/flux_shifts_OffAxis2023.root)) 
containing re-evaluated focussing, alignment and target condition systematics. 
The auxilliary script [`combine.C`](combine.C) contains the logic for creating 
our single input file. Note that the hadron production systematics were stored 
in a non-standard ROOT TH2 format and require extra code to extract, 
`combine.C` converts them to standard root objects in an equivalent format to 
the focussing, alignment, and target condition inputs.

## Usage in Samples

The main interface to the systematic parameter responses is through the 
[`OffAxisFluxUncertaintyHelper`](OffAxisFluxUncertaintyHelper.h) class.

*At initialisation time*, when events are read in, samples should retrieve two 
systematic bin identifiers for each event, one for the focussing, alignment, 
and target condition parameters and one for the hadron production parameters. 
A 'neutrino configuration' integer is used to identify the beam mode, neutrino 
species, and detector and should also be stored with the event to avoid 
re-calculation at step time. An example snippet for doing so is shown below:

```c++
//for each event, i
 duneobj->flux_syst_nu_config.push_back(
        flux_helper->GetNuConfig(duneobj->nupdgUnosc[i], isND, isFHC));

    // xdir in detsim and flux syst inputs are opposite.
    double syst_xpos_m = -(_det_x + _vtx_x)/100.0;
    duneobj->flux_focussing_syst_bin.push_back(flux_helper->GetFocussingBin(
        duneobj->flux_syst_nu_config.back(), enu_GeV, syst_xpos_m));
    duneobj->flux_hadprod_syst_bin.push_back(flux_helper->GetHadProdBin(
        duneobj->flux_syst_nu_config.back(), enu_GeV, syst_xpos_m));
```

The exact names of various variables may change in your particular sample, 
but hopefully, the usage is clear.

*At step time*, when weights should be retrieved, the previously determined 
neutrino configuration and systematic bin can be proferred to the interface, 
along with parameter identifiers and values, in exchange for response weights 
to be used in predicting varied event rates. An example snippet is shown 
below:

```c++
for (int i = 0; i < int(flux_helper->GetNFocussingParams()); i++) {
  //...
  double w = flux_helper->GetFluxFocussingWeight(
                i, *par, dunemcSamples[iSample].flux_syst_nu_config[iEvent],
                dunemcSamples[iSample].flux_focussing_syst_bin[iEvent]);
  //...
}
//...
for (int i = 0; i < int(flux_helper->GetNHadProdPCAComponents()); i++) {
  //...
  double w = flux_helper->GetFluxHadProdWeight(
                i, *par, dunemcSamples[iSample].flux_syst_nu_config[iEvent],
                dunemcSamples[iSample].flux_hadprod_syst_bin[iEvent]);
  //...
}
```

The response weight from each parameter represents independent variations and 
so can safetly multiplied together to produce a total event weight for some 
multi-parameter variation.

## Systematic Parameter Configuration

A `yaml` file containing the parameter configurations can be found in 
[FluxParameters_FD_and_PRISM.yaml](FluxParameters_FD_and_PRISM.yaml). A python 
script [`emit_flux_yaml.py`](emit_flux_yaml.py) is provided to emit this 
configuration given the input file of varied flux ratios 
([flux_variations_FD_and_PRISM_2023.root](flux_variations_FD_and_PRISM_2023.root)).

## Validations

The [`OffAxisFluxUncertaintyHelper`](OffAxisFluxUncertaintyHelper.h) class has 
been validated in MaCh3_DUNE by comparing the histograms directly from the 
input file with ratios built via direct calling of the interface code (outside 
of MaCh3) with the results of running [SigmaVariation](SigmaVariation.cpp) 
inside of MaCh3. At the time of committing, these validations passed for ND on 
axis and one off axis position. A script is provided to emit validation plots w
hen proffered an output file from a SigmaVariation run and the output of 
running the interface testing macro [`test_read_flux_systs.C`](test_read_flux_systs.C).

A pdf containing validation plots can be found here: 
[FluxVariationsValidation](FluxVariationsValidation.pdf).

## Auxilliary Scripts

* [`combine.C`](combine.C)
  + A script to combine two sets of input files from CAFAna into a single set 
  of inputs used by this library. Requires 
  [TH2Jagged](https://github.com/luketpickering/TH2Jagged) to be built as a 
  subdirectory of the directory that it is run in.
* [`test_read_flux_systs.C`](test_read_flux_systs.C)
  + A script to test that the interface code can read the input histograms. 
  Builds the interface code standalone with cling, and then loops over all the 
  parameters writing out some binning and parameter name information to prove 
  that it's managed to read the inputs. Also outputs a root file containing some 
  flux ratios retrieved via the interface to test that the right inputs are 
  being read for the right physics quantities. Used by `plot_flux_sigvar.py`.
* [`emit_flux_yaml.py`](emit_flux_yaml.py)
  + A script to read the file containing the input histograms and dump out a 
  systematic parameter `yaml` file.
* [`plot_flux_sigvar.py`](plot_flux_sigvar.py)
  + Validation plotting script, described in [Validations](#validations).
