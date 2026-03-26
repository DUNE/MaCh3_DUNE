#include "Samples/SampleHandlerBeamOffAxis.h"

#include "Samples/BeamOffAxis/ReadEvents.h"
#include "Samples/BeamOffAxis/Systematics.h"

#include <iostream>

#include <fstream>

namespace dune::beamoffaxis {

SampleHandlerBeamOffAxis::SampleHandlerBeamOffAxis(
    std::string mc_version_, ParameterHandlerGeneric *ParHandler_,
    const std::shared_ptr<OscillationHandler> &Oscillator_)
    : SampleHandlerFD(mc_version_, ParHandler_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();
}

void SampleHandlerBeamOffAxis::Init() {
  subsample_analysispot.resize(GetNsamples());
  subsample_is_numode.resize(GetNsamples());

  for (int iSubSample = 0; iSubSample < GetNsamples(); iSubSample++) {
    auto const &sample_conf = SampleManager->raw()[GetSampleTitle(iSubSample)];
    subsample_analysispot[iSubSample] =
        Get<double>(sample_conf["POT"], __FILE__, __LINE__);
    subsample_is_numode[iSubSample] =
        Get<bool>(sample_conf["is_numode"], __FILE__, __LINE__);
  }
}

void SampleHandlerBeamOffAxis::SetupSplines() {

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

void SampleHandlerBeamOffAxis::RegisterFunctionalParameters() {

  RegisterIndividualFunctionalParameter(DUNEMCEvents,
                                        {
                                            "ContainedMuonEnergyScale",
                                            "ContainedMuonSqrtEnergyScale",
                                            "ContainedMuonInvSqrtEnergyScale",
                                            "TrackedMuonEnergyScale",
                                            "TrackedMuonSqrtEnergyScale",
                                            "TrackedMuonInvSqrtEnergyScale",
                                            "EMEnergyScale",
                                            "EMSqrtEnergyScale",
                                            "EMInvSqrtEnergyScale",
                                            "ChgHadEnergyScale",
                                            "ChgHadSqrtEnergyScale",
                                            "ChgHadInvSqrtEnergyScale",
                                            "NeutronEnergyScale",
                                            "NeutronSqrtEnergyScale",
                                            "NeutronInvSqrtEnergyScale",
                                            "TotalEnergyScale",
                                            "TotalSqrtEnergyScale",
                                            "TotalInvSqrtEnergyScale",
                                        },
                                        EnergyScales);

  RegisterIndividualFunctionalParameter(
      DUNEMCEvents,
      {"MuonEnergyResolution", "EMEnergyResolution", "ChgHadEnergyResolution",
       "NeutronEnergyResolution"},
      ParticleEnergyResolutions);

  if (ParHandler->GetNumParFromGroup("Flux")) {
    RegisterIndividualFunctionalParameter(
        DUNEMCEvents, GetFluxFocussingParamNames(), UpdateFluxFocussingWeight);

    RegisterIndividualFunctionalParameter(
        DUNEMCEvents, GetFluxHadProdParamNames(), UpdateFluxHadProdWeight);
  }
}

void SampleHandlerBeamOffAxis::ResetShifts(int iEvent) {
  DUNEMCEvents[iEvent].varied_reco = DUNEMCEvents[iEvent].reco;
  DUNEMCEvents[iEvent].syst.flux.total_weight = 1.0;
}

void SampleHandlerBeamOffAxis::FinaliseShifts(int iEvent) {
  CalculateVariedCompositeQuantities(DUNEMCEvents[iEvent]);
}

void SampleHandlerBeamOffAxis::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < DUNEMCEvents.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(
        &(DUNEMCEvents[i].weights.pot));
    MCEvents[i].total_weight_pointers.push_back(
        &(DUNEMCEvents[i].syst.flux.total_weight));
  }
}

int SampleHandlerBeamOffAxis::SetupExperimentMC() {

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");

  bool do_flux_systematics = ParHandler->GetNumParFromGroup("Flux");

  for (int iSubSample = 0; iSubSample < int(SampleDetails.size());
       iSubSample++) {
    MACH3LOG_INFO("-- subsample[{}]: {} ", iSubSample,
                  SampleDetails[iSubSample].SampleTitle);

    TChain MetaChain("meta");
    TChain CAFChain("cafTree");
    for (const std::string &filename : SampleDetails[iSubSample].mc_files) {
      MACH3LOG_INFO("-- -- Adding file to TChain: {}", filename);
      if (!CAFChain.Add(filename.c_str(), -1)) {
        MACH3LOG_ERROR("Could not add file {} to TChain, please check the file "
                       "exists and is readable",
                       filename);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      MetaChain.Add(filename.c_str(), -1);
    }

    double subsample_cafpot = GetPOT(MetaChain);
    auto sample_evs = ReadEvents(CAFChain);

    // fix up any analysis specific information
    for (auto &ev : sample_evs) {

      ev.subsample = iSubSample;
      ev.is_numode = subsample_is_numode[iSubSample];

      ev.weights.pot = subsample_analysispot[iSubSample] / subsample_cafpot;

      ev.truth.mach3_mode =
          Modes->GetModeFromGenerator(std::abs(ev.truth.generator_mode));
      if (!ev.truth.is_cc) {
        // Account for no ability to distinguish CC/NC
        ev.truth.mach3_mode += 14;
      }
      if (ev.truth.mach3_mode > 15) {
        // Account for no NCSingleKaon
        ev.truth.mach3_mode -= 1;
      }

      ev.syst.flux.total_weight = 1;
      if (do_flux_systematics) {
        std::tie(ev.syst.flux.focussing_ratio, ev.syst.flux.hadprod_ratio) =
            GetFluxVariationRatios(ev.truth.nu.pdg, ev.truth.nu.e,
                                   ev.truth.vtx.off_axis_pos_m, true);
      }
    }

    DUNEMCEvents.reserve(DUNEMCEvents.size() + CAFChain.GetEntries());
    std::copy(sample_evs.begin(), sample_evs.end(),
              std::back_inserter(DUNEMCEvents));
  }

  return int(DUNEMCEvents.size());
}

const double *
SampleHandlerBeamOffAxis::GetPointerToKinematicParameter(KinematicTypes KinPar,
                                                         int iEvent) {
  return ResolveKinematicEventMember(KinPar, DUNEMCEvents[iEvent]);
}

const double *SampleHandlerBeamOffAxis::GetPointerToKinematicParameter(
    std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(
      ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

const double *SampleHandlerBeamOffAxis::GetPointerToKinematicParameter(
    double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamOffAxis::ReturnKinematicParameter(int KinematicVariable,
                                                          int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamOffAxis::ReturnKinematicParameter(
    std::string KinematicParameter, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iEvent);
}

void SampleHandlerBeamOffAxis::SetupFDMC() {
  size_t iEvent = 0;
  for (auto const &ev : DUNEMCEvents) {
    MCEvents[iEvent].Target = ev.truth.target_a;
    MCEvents[iEvent].mode = ev.truth.mach3_mode;
    MCEvents[iEvent].isNC = !ev.truth.is_cc;

    MCEvents[iEvent].enu_true = ev.truth.nu.e;
    MCEvents[iEvent].nupdg = ev.truth.nu.pdg;
    MCEvents[iEvent].nupdgUnosc = ev.truth.nu.pdg_unosc;

    MCEvents[iEvent].NominalSample = ev.subsample;

    iEvent++;
  }
}

} // namespace dune::beamoffaxis
