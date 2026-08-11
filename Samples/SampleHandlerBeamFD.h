#ifndef _samplePDFDUNEBeamFD_h_
#define _samplePDFDUNEBeamFD_h_

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerBase.h"

#include "StructsDUNE.h"
/// @brief Base class for handling FD Beam samples
class SampleHandlerBeamFD : virtual public SampleHandlerBase
{
public:


  /// @brief SampleHandler FD beam Constructor
  /// @param mc_version Config Name
  /// @param xsec_cov Cross-section covariance matrix
  /// @param osc_cov Oscillation covariance matrix
  /// @param Oscillator_ Shared Oscillation Handler object
  SampleHandlerBeamFD(std::string mc_version, ParameterHandlerGeneric* xsec_cov, const std::shared_ptr<OscillationHandler>&  Oscillator_);

  /// @brief destructor
  ~SampleHandlerBeamFD();

  /// @brief Enum to identify kinematics
  enum KinematicTypes
  {
    kTrueNeutrinoEnergy,
    kRecoNeutrinoEnergy,
    kTrueXPos,
    kTrueYPos,
    kTrueZPos,
    kCVNNumu,
    kCVNNue,
    kM3Mode,
    kOscChannel,
    kIsFHC,
    kTargetNucleus,
    kTrueCCnue,
    kTrueCCnumu
  };

protected:
  /// @brief Initialises object
  void Init();

  /// @brief Initialise data hist (can be overridden)
  void InititialiseData();

  /// @brief Function to setup MC from file
  /// @return Total number of events
  int SetupExperimentMC();

  /// @brief Tells FD base which variables to point to/be set to
  void SetupMC();

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void AddAdditionalWeightPointers();
  void SetupSplines();

  void RegisterFunctionalParameters() override;
  void ResetShifts(int iEvent) override;

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as int
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter (const int KinematicVariable, const int iEvent) const override;

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter as double (gets cast -> int)
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(const int KinematicVariable, const int iEvent) const override;

  //DB functions which could be initialised to do something which is non-trivial
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iEvent) {(void)iEvent; return 1.;}

  // dunemc
  /// DUNE MC samples
  std::vector<dunemc_beamfd> dunemcSamples;
  std::vector<BeamFDSampleInfo> beamFDSampleDetails;

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
    {"TrueXPos",kTrueXPos},
    {"TrueYPos",kTrueYPos},
    {"TrueZPos",kTrueZPos},
    {"CVNNumu",kCVNNumu},
    {"CVNNue",kCVNNue},
    {"Mode",kM3Mode},
    {"OscillationChannel",kOscChannel},
    {"IsFHC",kIsFHC},
    {"TargetNucleus", kTargetNucleus},
    {"IsTrueCCnue", kTrueCCnue},
    {"IsTrueCCnumu", kTrueCCnumu}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kTrueXPos,"TrueXPos"},
    {kTrueYPos,"TrueYPos"},
    {kTrueZPos,"TrueZPos"},
    {kCVNNumu,"CVNNumu"},
    {kCVNNue,"CVNNue"},
    {kM3Mode,"Mode"},
    {kOscChannel,"OscillationChannel"},
    {kIsFHC,"IsFHC"},
    {kTargetNucleus, "TargetNucleus"},
    {kTrueCCnue,"IsTrueCCnue"},
    {kTrueCCnumu,"IsTrueCCnumu"}
  };
  std::unordered_map<std::string, std::vector<double>> norm_map;

  /// @brief Downsampling step (e.g. 10 means only every 10th event is used in the fit). Default is 1 (no downsampling).
  unsigned int downsamplingStep;

  /// @brief Cleanup memory
  void CleanMemoryBeforeFit() override {};
};



#endif
