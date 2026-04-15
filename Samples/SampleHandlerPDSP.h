#pragma once

#include "Samples/SampleHandlerFD.h"
#include "StructsPDSP.h"
#include <random>

class SampleHandlerPDSP : public SampleHandlerFD
{
 public:
  SampleHandlerPDSP(const std::string& config_name, ParameterHandlerGeneric* parameter_handler);
  virtual ~SampleHandlerPDSP();

  enum KinematicTypes {kTrueKEInt, kRecoKEInt};
  
  // =============================================

 protected:
  void Init() override;

  ///@brief Setup our spline file, this calls InitialseSplineObject() under the hood
  void SetupSplines() override;

  int SetupExperimentMC() override;

  void CleanMemoryBeforeFit() override;


  double ReturnKinematicParameter(KinematicTypes KinPar, int iEvent);
  double ReturnKinematicParameter(int KinematicVariable, int iEvent) override;
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent) override;
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) override;
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) override;

  void SetupFDMC() override;
  void CalcWeightFunc(int iEvent) override {return; (void)iEvent;}

  std::vector<PDSPMCInfo> PDSPSamples;
  std::vector<PDSPMCPlottingInfo> PDSPPlottingSamples;

  const std::unordered_map<std::string, int> KinematicParametersPDSP = {
    {"TrueKEInt", kTrueKEInt},
    {"RecoKEInt", kRecoKEInt},
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersPDSP = {
    {kTrueKEInt, "TrueKEInt"},
    {kRecoKEInt, "RecoKEInt"},
  };

  // functional parameters, currently have none for the time being
  void RegisterFunctionalParameters() override;
};
