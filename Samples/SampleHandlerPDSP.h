#pragma once

#include "Samples/SampleHandlerBase.h"
#include "StructsPDSP.h"
#include <random>

class SampleHandlerPDSP : virtual public SampleHandlerBase
{
 public:
  SampleHandlerPDSP(const std::string& config_name, ParameterHandlerGeneric* parameter_handler);
  virtual ~SampleHandlerPDSP();

  enum KinematicTypes {kTrueKEIni, kTrueKEInt, kRecoKEIni, kRecoKEInt, kMode, kOscChannel, kTargetNucleus, kTrueEndZ, kRecoEndZ};

  TH1* GetDataHistogramFromInputs(const int Sample) const;
  
  // =============================================

 protected:
  void Init() override;

  ///@brief Setup our spline file, this calls InitialseSplineObject() under the hood
  void SetupSplines() override;

  int SetupExperimentMC() override;

  void SetupMC() override;

  void InititialiseData() override;

  void CleanMemoryBeforeFit() override;

  void AddAdditionalWeightPointers() override;

  double ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const override;
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent) const;
  const double* GetPointerToKinematicParameter(const int KinematicVariable, const int iEvent) const override;

  void CalcWeightFunc(const int iEvent) override {return; (void)iEvent;}

  std::vector<MetaData> PDSPSampleMetaData;
  std::vector<PDSPMCInfo> PDSPSamples;
  std::vector<PDSPMCPlottingInfo> PDSPPlottingSamples;

  const std::unordered_map<std::string, int> KinematicParametersPDSP = {
    {"TrueKEIni", kTrueKEIni},
    {"TrueKEInt", kTrueKEInt},
    {"RecoKEIni", kRecoKEIni},
    {"RecoKEInt", kRecoKEInt},
    {"Mode", kMode},
    {"OscillationChannel", kOscChannel},
    {"TargetNucleus", kTargetNucleus},
    {"TrueEndZ", kTrueEndZ},
    {"RecoEndZ", kRecoEndZ}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersPDSP = {
    {kTrueKEIni, "TrueKEIni"},
    {kTrueKEInt, "TrueKEInt"},
    {kRecoKEIni, "RecoKEIni"},
    {kRecoKEInt, "RecoKEInt"},
    {kMode, "Mode"},
    {kOscChannel, "OscillationChannel"},
    {kTargetNucleus, "TargetNucleus"},
    {kTrueEndZ, "TrueEndZ"},
    {kRecoEndZ, "RecoEndZ"}
  };

  // functional parameters, currently have none for the time being
  void RegisterFunctionalParameters() override;

  M3::float_t MCGlobalScale;
};
