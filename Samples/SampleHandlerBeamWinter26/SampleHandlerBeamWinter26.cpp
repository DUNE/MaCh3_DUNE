#include "Samples/SampleHandlerBeamWinter26/SampleHandlerBeamWinter26.h"

#include "Samples/SampleHandlerBeamWinter26/ReadEvents.h"

#include "Splines/BinnedSplineHandlerDUNE.h"

#include <iostream>

#include <fstream>

namespace dune::beam {

template <bool READ_EVENT_SPLINES>
SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::SampleHandlerBeamWinter26(
    std::string mc_version_, ParameterHandlerGeneric *ParHandler_,
    const std::shared_ptr<OscillationHandler> &Oscillator_)
    : SampleHandlerBase(mc_version_, ParHandler_, Oscillator_) {

  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();
}

template <bool READ_EVENT_SPLINES>
void SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::Init() {}

template <bool READ_EVENT_SPLINES>
void SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::SetupSplines() {
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

template <bool READ_EVENT_SPLINES>
void SampleHandlerBeamWinter26<
    READ_EVENT_SPLINES>::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < DUNEMCEvents.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(DUNEMCEvents[i].weights.pot));
  }
}

template <bool READ_EVENT_SPLINES>
int SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::SetupExperimentMC() {

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");

  std::vector<std::array<float, 2>> RecoSampleRanges;

  auto sample_name = Get<std::vector<std::string>>(
      SampleManager->raw()["Samples"], __FILE__, __LINE__);

  for (int i = 0; i < GetNSamples(); i++) {
    auto first_cut = Get<YAML::Node>(
        SampleManager->raw()[sample_name[i]]["SelectionCuts"][0], __FILE__,
        __LINE__);
    auto kinstr =
        Get<std::string>(first_cut["KinematicStr"], __FILE__, __LINE__);

    if (kinstr == ReversedKinematicParametersDUNE.at(kRecoSample)) {
      RecoSampleRanges.push_back(
          Get<std::array<float, 2>>(first_cut["Bounds"], __FILE__, __LINE__));
    } else {
      MACH3LOG_ERROR("Expected to only find a single Selection cut cutting on "
                     "{}. But found a cut on {}",
                     ReversedKinematicParametersDUNE.at(kRecoSample), kinstr);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  std::map<std::string, std::vector<std::pair<std::string, float>>>
      input_mc_event_descriptors;

  for (auto const &file_descriptor :
       SampleManager->raw()["InputFiles"]["MCEvents"]) {

    auto tag = Get<std::string>(file_descriptor["Tag"], __FILE__, __LINE__);
    auto file_location =
        Get<std::string>(file_descriptor["FileLocation"], __FILE__, __LINE__);
    auto downsamplefraction = GetFromManager<float>(
        file_descriptor["DownsampleFraction"], 0.0, __FILE__, __LINE__);

    input_mc_event_descriptors[tag].push_back(
        std::make_pair(file_location, downsamplefraction));

    MACH3LOG_INFO("-- Found input event descriptor: Tag: {}, FileLocation: {}, "
                  "DownsampleFraction: {}",
                  tag, file_location, downsamplefraction);
  }

  for (auto const &[tag, input_files] : input_mc_event_descriptors) {

    float tag_pot =
        Get<float>(SampleManager->raw()["POT"][tag], __FILE__, __LINE__);
    size_t tag_id = mc_tags.size();
    mc_tags.push_back(MCTag{tag, tag_pot});

    bool is_numode = tag.find("numode") != std::string::npos;

    float tag_input_pot = 0;

    for (auto const &[filename, downsamplefraction] : input_files) {
      if (filename.empty()) {
        MACH3LOG_INFO("-- -- Skipping empty filename entry");
        continue;
      }
      std::unique_ptr<TFile> metaf(TFile::Open(filename.c_str(), "READ"));
      float file_descriptor_pot = GetPOT(*metaf);
      tag_input_pot += file_descriptor_pot * (1 - downsamplefraction);
      MACH3LOG_INFO(
          "-- -- Read {:.3G} input POT (with downsample weight of: {:.2f})",
          file_descriptor_pot, (1 - downsamplefraction));
    }

    MACH3LOG_INFO("-- Read {:.3G} total POT for tag: {}, which has analysis "
                  "POT of {:.3G}",
                  tag_input_pot, tag, tag_pot);

    for (auto const &[filename, downsamplefraction] : input_files) {
      if (filename.empty()) {
        MACH3LOG_INFO("-- -- Skipping empty filename entry");
        continue;
      }

      std::unique_ptr<TFile> eventf(TFile::Open(filename.c_str(), "READ"));

      auto [sample_evs, entries] =
          ReadEvents<READ_EVENT_SPLINES>(*eventf, downsamplefraction);
      MACH3LOG_INFO(
          "-- -- Read: {}/{} events (with downsample fraction of: {:.2f})",
          sample_evs.size(), entries, downsamplefraction);

      // fix up any analysis specific information
      for (auto &ev : sample_evs) {

        ev.tag_id = tag_id;
        ev.is_numode = is_numode;

        // here we manually apply the specified sample cut to ignore events that
        // didn't enter a selection.
        ev.sample = -1;
        for (size_t i = 0; i < RecoSampleRanges.size(); ++i) {
          if ((ev.reco.sample > RecoSampleRanges[i][0]) &&
              (ev.reco.sample < RecoSampleRanges[i][1])) {
            ev.sample = int(i);
          }
        }

        if (ev.sample < 0) {
          continue;
        }

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

        ev.weights.pot = tag_pot / tag_input_pot;

        DUNEMCEvents.emplace_back(std::move(ev));
      }
    }
  }

  return int(DUNEMCEvents.size());
}

template <bool READ_EVENT_SPLINES>
void SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::SetupMC() {

  size_t iEvent = 0;
  for (auto const &ev : DUNEMCEvents) {
    MCEvents[iEvent].isNC = !ev.truth.is_cc;

    MCEvents[iEvent].enu_true = ev.truth.nu.e;
    MCEvents[iEvent].nupdg = ev.truth.nu.pdg;
    MCEvents[iEvent].nupdgUnosc = ev.truth.nu.pdg_unosc;

    MCEvents[iEvent].NominalSample = ev.sample;

    iEvent++;
  }
}

template <bool READ_EVENT_SPLINES>
void SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::InititialiseData() {
  // Reweight MC to match
  Reweight();
  // set asimov data
  for (int iSample = 0; iSample < GetNSamples(); iSample++) {
    AddData(iSample, GetMCArray(iSample));
  }
}

template <bool READ_EVENT_SPLINES>
const double *
SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::GetPointerToKinematicParameter(
    int KinematicVariable, int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return ResolveKinematicEventMember(KinPar, DUNEMCEvents[iEvent]);
}

template <bool READ_EVENT_SPLINES>
double SampleHandlerBeamWinter26<READ_EVENT_SPLINES>::ReturnKinematicParameter(
    int KinematicVariable, int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

template class SampleHandlerBeamWinter26<true>;
template class SampleHandlerBeamWinter26<false>;

} // namespace dune::beam
