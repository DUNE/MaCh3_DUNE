#pragma once

#include "Samples/SampleHandlerBeamWinter26/EventInfo.h"
#include "Samples/SampleHandlerBeamWinter26/Projections.h"

#include "Samples/SampleHandlerBase.h"

namespace dune::beam {

template <bool READ_EVENT_SPLINES = false>
class SampleHandlerBeamWinter26 : public SampleHandlerBase {

public:
  SampleHandlerBeamWinter26(
      std::string mc_version, ParameterHandlerGeneric *xsec_cov,
      const std::shared_ptr<OscillationHandler> &Oscillator);

  ~SampleHandlerBeamWinter26() {}

  using TEvent =
      std::conditional_t<READ_EVENT_SPLINES, EventInfoSplines, EventInfo>;

  std::vector<std::pair<std::string, std::vector<float>>> EventSplineMetadata;
  std::vector<TEvent> DUNEMCEvents;

  struct MCTag {
    std::string name;
    float POT;
  };
  std::vector<MCTag> mc_tags;

  int GetTagIndex(std::string const &tag) {
    int index = 0;
    for (auto const &t : mc_tags) {
      if (t.name == tag) {
        return index;
      }
      index++;
    }
    return -1;
  }

  void CleanMemoryBeforeFit() {}

  // if READ_EVENT_SPLINES == true, apply knot response kt_i for param pr_i
  void SetEventResponseWeights(int pr_i = -1, int kt_i = -1);

  auto const &GetBinning() { return Binning; }

protected:
  void Init() override;
  int SetupExperimentMC() override;
  void SetupMC() override;
  void AddAdditionalWeightPointers() override;
  void SetupSplines() override;
  void InititialiseData();

  double ReturnKinematicParameter(int KinematicVariable, int iEvent) const;

  const double *GetPointerToKinematicParameter(int KinematicVariable,
                                               int iEvent) const;
};
} // namespace dune::beam
