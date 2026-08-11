#pragma once

#include "Samples/BeamFDStandardRecord/EventInfo.h"
#include "Samples/BeamFDStandardRecord/Projections.h"
#include "Samples/SampleHandlerBase.h"

namespace dune::beamfd {

class SampleHandlerBeamFDStandardRecord : public SampleHandlerBase {

public:
  SampleHandlerBeamFDStandardRecord(
      std::string mc_version, ParameterHandlerGeneric *xsec_cov,
      const std::shared_ptr<OscillationHandler> &Oscillator);

  ~SampleHandlerBeamFDStandardRecord() {}

  std::vector<dune::beamfd::EventInfo> DUNEMCEvents;

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
  void RegisterFunctionalParameters() override;
  void ResetShifts(int iEvent) override;
  void InititialiseData();

  double ReturnKinematicParameter(int KinematicVariable, int iEvent) const;

  const double *GetPointerToKinematicParameter(int KinematicVariable,
                                               int iEvent) const;
};
} // namespace dune::beamfd
