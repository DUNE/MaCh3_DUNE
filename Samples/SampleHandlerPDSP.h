#pragma once

#include "Samples/SampleHandlerFD.h"
#include "StructsPDSP.h"
#include <random>

class SampleHandlerPDSP : public SampleHandlerFD
{
 public:
  SampleHandlerPDSP(const std::string& config_name, ParameterHandlerGeneric* parameter_handler);
  virtual ~SampleHandlerPDSP();

  enum KinematicTypes {kTrueKEInt, kRecoKEInt, kMode, kOscChannel, kRecoEndZ};
  
  // =============================================

 protected:
  void Init() override;

  ///@brief Setup our spline file, this calls InitialseSplineObject() under the hood
  void SetupSplines() override;

  int SetupExperimentMC() override;

  void CleanMemoryBeforeFit() override;

  void AddAdditionalWeightPointers() override;

  double ReturnKinematicParameter(KinematicTypes KinPar, int iEvent);
  double ReturnKinematicParameter(int KinematicVariable, int iEvent) override;
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent) override;
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) override;
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) override;

  void SetupFDMC() override;
  void CalcWeightFunc(int iEvent) override {return; (void)iEvent;}

  std::vector<MetaData> PDSPSampleMetaData;
  std::vector<PDSPMCInfo> PDSPSamples;
  std::vector<PDSPMCPlottingInfo> PDSPPlottingSamples;

  const std::unordered_map<std::string, int> KinematicParametersPDSP = {
    {"TrueKEInt", kTrueKEInt},
    {"RecoKEInt", kRecoKEInt},
    {"Mode", kMode},
    {"OscillationChannel", kOscChannel},
    {"RecoEndZ", kRecoEndZ}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersPDSP = {
    {kTrueKEInt, "TrueKEInt"},
    {kRecoKEInt, "RecoKEInt"},
    {kMode, "Mode"},
    {kOscChannel, "OscillationChannel"},
    {kRecoEndZ, "RecoEndZ"},
  };

  // functional parameters, currently have none for the time being
  void RegisterFunctionalParameters() override;

  // Placeholder for nupdg/nupdgUnosc/Target pointers — unused in PDSP but must not be null
  static const int DummyInt = 0;
};
