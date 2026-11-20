#pragma once

#include "Samples/SampleHandlerBeamWinter26/EventInfo.h"
#include "Samples/SampleHandlerBeamWinter26/Projections.h"

#include "Samples/SampleHandlerBase.h"

namespace dune::beam {

class SampleHandlerBeamWinter26 : public SampleHandlerBase {

public:
  SampleHandlerBeamWinter26(
      std::string mc_version, ParameterHandlerGeneric *xsec_cov,
      const std::shared_ptr<OscillationHandler> &Oscillator);

  ~SampleHandlerBeamWinter26() {}

  std::vector<dune::beam::EventInfo> DUNEMCEvents;

  struct MCTag {
    std::string name;
    float POT;
  };
  std::vector<MCTag> mc_tags;

  int GetTagIndex(std::string const &tag){
    int index = 0;
    for(auto const &t : mc_tags){
      if(t.name == tag){
        return index;
      }
      index++;
    }
    return -1;
  }

  void CleanMemoryBeforeFit() {}

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
