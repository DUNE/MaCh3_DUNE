#ifndef _SampleHandlerAtm_h_
#define _SampleHandlerAtm_h_

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerBase.h"

#include "StructsDUNE.h"
/// @brief Base class for handling atmospheric samples
class SampleHandlerAtm : virtual public SampleHandlerBase
{
public:
  /// @brief Constructor
  /// @param mc_version Configuration file
  /// @param xsec_cov cross-section covariance matrix
  /// @param osc_cov oscillation covariance matrix
  SampleHandlerAtm(std::string mc_version, ParameterHandlerGeneric* xsec_cov, const std::shared_ptr<OscillationHandler>&  Oscillator_);
  /// @brief destructor
  ~SampleHandlerAtm();

  /// @brief Enum to identify kinematics
  enum KinematicTypes
  {
    kTrueNeutrinoEnergy,
    kRecoNeutrinoEnergy,
    kTrueCosZ,
    kRecoCosZ,
    kOscChannel,
    kMode,
    kTargetNucleus
  };

protected:
  /// @brief Initialises object
  void Init();

  /// @brief Initialise data hist (can be overridden)
  void InititialiseData() override;

  /// @brief Function to setup MC from file
  /// @param iSample sample ID
  /// @return Total number of events
  int SetupExperimentMC();

  /// @brief Tells FD base which variables to point to/be set to
  /// @param iSample Sample ID
  void SetupMC();

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void AddAdditionalWeightPointers();

  /// @brief Sets up splines 
  void SetupSplines();

  /// @brief Cleanup memory
  void CleanMemoryBeforeFit() override {};

  void RegisterFunctionalParameters() override {};
  
  //DB functions which could be initialised to do something which is non-trivial
  
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iSample sample ID
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iEvent) {(void) iEvent; return 1.;}
  
  /// @brief NOT IMPLEMENTED: Apply kinematic shifts
  /// @param iSample Sample Number
  /// @param iEvent Event number
  void applyShifts(int iEvent) {(void)iEvent;}
  
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(const int KinematicVariable, const int iEvent) const override;

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const override;

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
    {"TrueCosineZ",kTrueCosZ},
    {"RecoCosineZ",kRecoCosZ},
    {"OscillationChannel",kOscChannel},
    {"Mode",kMode},
    {"TargetNucleus", kTargetNucleus}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kTrueCosZ,"TrueCosineZ"},    
    {kRecoCosZ,"RecoCosineZ"},
    {kOscChannel,"OscillationChannel"},
    {kMode,"Mode"},
    {kTargetNucleus, "TargetNucleus"},
  };
  
  /// Array filled with MC samples for each oscillation channel
  std::vector<dunemc_atm> dunemcSamples;

  /// Is the sample e-like
  std::vector<int> IsELike;

  /// Multiplicative scaling to scale from the assumed 400ktyr value in the CAF files
  double ExposureScaling;
};

#endif
