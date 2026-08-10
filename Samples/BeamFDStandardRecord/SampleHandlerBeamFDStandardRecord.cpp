#include "Samples/BeamFDStandardRecord/SampleHandlerBeamFDStandardRecord.h"

#include "Samples/BeamFDStandardRecord/ReadEvents.h"

#include "Systematics/Beam.h"

#include "Splines/BinnedSplineHandlerDUNE.h"

#include <iostream>

#include <fstream>

namespace dune::beamfd {

SampleHandlerBeamFDStandardRecord::SampleHandlerBeamFDStandardRecord(
    std::string mc_version_, ParameterHandlerGeneric *ParHandler_,
    const std::shared_ptr<OscillationHandler> &Oscillator_)
    : SampleHandlerBase(mc_version_, ParHandler_, Oscillator_) {

  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();
}

void SampleHandlerBeamFDStandardRecord::Init() {
  subsample_analysispot.resize(GetNSamples());
  subsample_is_numode.resize(GetNSamples());

  for (int iSubSample = 0; iSubSample < GetNSamples(); iSubSample++) {
    auto const &sample_conf = SampleManager->raw()[GetSampleTitle(iSubSample)];
    subsample_analysispot[iSubSample] =
        Get<double>(sample_conf["POT"], __FILE__, __LINE__);
    subsample_is_numode[iSubSample] =
        Get<bool>(sample_conf["is_numode"], __FILE__, __LINE__);
  }
}

void SampleHandlerBeamFDStandardRecord::SetupSplines() {
  if (!ParHandler) {
    return;
  }

  ///@todo move all of the spline setup into core
  int num_splines = 0;
  for (int iSubSample = 0; iSubSample < int(SampleDetails.size());
       iSubSample++) {
    num_splines += ParHandler->GetNumParamsFromSampleName(
        GetSampleTitle(iSubSample), kSpline);
  }

  if (num_splines > 0) {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will create a spline object",
        num_splines);
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(
        new BinnedSplineHandlerDUNE(ParHandler, Modes.get()));
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will not load or "
        "evaluate splines",
        ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = nullptr;
  }
}

void SampleHandlerBeamFDStandardRecord::RegisterFunctionalParameters() {
  if (!ParHandler) {
    return;
  }

  if (ParHandler->GetNumParFromGroup("Flux")) {
    RegisterIndividualFunctionalParameter(
        DUNEMCEvents, syst::GetFluxFocussingParamNames(),
        [](std::vector<double> const &par_vals, EventInfo &ev) {
          for (size_t i = 0; i < par_vals.size(); ++i) {
            ev.syst.flux.total_weight *=
                1 + (par_vals[i] * ev.syst.flux.focussing_ratio[i]);
          }
        });

    RegisterIndividualFunctionalParameter(
        DUNEMCEvents, syst::GetFluxHadProdParamNames(),
        [](std::vector<double> const &par_vals, EventInfo &ev) {
          for (size_t i = 0; i < par_vals.size(); ++i) {
            ev.syst.flux.total_weight *=
                1 + (par_vals[i] * ev.syst.flux.hadprod_ratio[i]);
          }
        });
  }
}

void SampleHandlerBeamFDStandardRecord::ResetShifts(int iEvent) {
  auto &ev = DUNEMCEvents[iEvent];

  // flux weights
  ev.syst.flux.total_weight = 1.0;
}

void SampleHandlerBeamFDStandardRecord::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < DUNEMCEvents.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(DUNEMCEvents[i].weights.pot));
    MCEvents[i].total_weight_pointers.push_back(
        &(DUNEMCEvents[i].syst.flux.total_weight));
  }
}

int SampleHandlerBeamFDStandardRecord::SetupExperimentMC() {

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");

  bool do_flux_systematics =
      ParHandler && ParHandler->GetNumParFromGroup("Flux");

  for (int iSubSample = 0; iSubSample < int(SampleDetails.size());
       iSubSample++) {
    MACH3LOG_INFO("-- subsample[{}]: {} ", iSubSample,
                  SampleDetails[iSubSample].SampleTitle);

    TChain MetaChain("meta");
    TChain CAFChain("cafTree");
    for (auto const &osc_channel_filenames :
         SampleDetails[iSubSample].mc_files) {
      for (const std::string &filename : osc_channel_filenames) {
        if (filename.empty()) {
          MACH3LOG_INFO("-- -- Skipping empty filename entry");
          continue;
        }
        MACH3LOG_INFO("-- -- Adding file to TChain: {}", filename);
        if (!CAFChain.Add(filename.c_str(), -1)) {
          MACH3LOG_ERROR(
              "Could not add file {} to TChain, please check the file "
              "exists and is readable",
              filename);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        MetaChain.Add(filename.c_str(), -1);
      }
    }

    double subsample_cafpot = GetPOT(MetaChain);
    auto sample_evs = ReadEvents(CAFChain);

    // fix up any analysis specific information
    for (auto &ev : sample_evs) {

      ev.subsample = iSubSample;
      ev.is_numode = subsample_is_numode[iSubSample];

      ev.weights.pot = subsample_analysispot[iSubSample] / subsample_cafpot;

      ev.syst.flux.total_weight = 1;
      if (do_flux_systematics) {
        // do stuff here
      }
    }

    DUNEMCEvents.reserve(DUNEMCEvents.size() + sample_evs.size());
    std::copy(sample_evs.begin(), sample_evs.end(),
              std::back_inserter(DUNEMCEvents));
  }

  return int(DUNEMCEvents.size());
}

void SampleHandlerBeamFDStandardRecord::SetupMC() {
  size_t iEvent = 0;
  for (auto const &ev : DUNEMCEvents) {
    MCEvents[iEvent].isNC = !ev.truth.is_cc;

    MCEvents[iEvent].enu_true = ev.truth.nu.e;
    MCEvents[iEvent].nupdg = ev.truth.nu.pdg;
    MCEvents[iEvent].nupdgUnosc = ev.truth.nu.pdg_unosc;

    MCEvents[iEvent].NominalSample = ev.subsample;

    iEvent++;
  }
}

void SampleHandlerBeamFDStandardRecord::InititialiseData() {
  // Reweight MC to match
  Reweight();
  // set asimov data
  for (int iSample = 0; iSample < GetNSamples(); iSample++) {
    AddData(iSample, GetMCArray(iSample));
  }
}

const double *SampleHandlerBeamFDStandardRecord::GetPointerToKinematicParameter(
    int KinematicVariable, int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamFDStandardRecord::ReturnKinematicParameter(
    int KinematicVariable, int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

} // namespace dune::beamfd
