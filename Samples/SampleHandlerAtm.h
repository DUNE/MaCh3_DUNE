#ifndef _SampleHandlerAtm_h_
#define _SampleHandlerAtm_h_

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerFD.h"

#include "StructsDUNE.h"
/// @brief Base class for handling atmospheric samples
class SampleHandlerAtm : virtual public SampleHandlerFD
{
public:
  /// @brief Constructor
  /// @param mc_version Configuration file
  /// @param xsec_cov cross-section covariance matrix
  /// @param osc_cov oscillation covariance matrix
  SampleHandlerAtm(std::string mc_version, ParameterHandlerGeneric* xsec_cov, ParameterHandlerOsc* osc_cov);
  /// @brief destructor
  ~SampleHandlerAtm();

  /// @brief Enum to identify kinematics
  enum KinematicTypes{kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueCosZ,kRecoCosZ,kOscChannel,kMode};
  
protected:
  /// @brief Initialises object
  void Init();

  /// @brief Function to setup MC from file
  /// @param iSample sample ID
  /// @return Total number of events
  int SetupExperimentMC(int iSample);

  /// @brief Tells FD base which variables to point to/be set to
  /// @param iSample Sample ID
  void SetupFDMC(int iSample);

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void SetupWeightPointers();

  /// @brief Sets up splines 
  void SetupSplines();

  void RegisterFunctionalParameters() override {};
  
  //DB functions which could be initialised to do something which is non-trivial
  
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iSample sample ID
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iSample, int iEvent) {(void)iSample; (void) iEvent; return 1.;}
  
  /// @brief NOT IMPLEMENTED: Apply kinematic shifts
  /// @param iSample Sample Number
  /// @param iEvent Event number
  void applyShifts(int iSample, int iEvent) {(void)iSample; (void)iEvent;}
  
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinPar Kinematic parameter enum val
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(int KinematicVariable, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  /// @brief Gets binning for a given parameter
  /// @param KinematicParameterStr Parameter name
  /// @return Vector containing parameter bins
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameterStr);

  /// @brief Gets binning for a given parameter
  /// @param KinPar Parameter ID
  /// @return Vector containing parameter bins
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinPar);
  
  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
    {"TrueCosineZ",kTrueCosZ},
    {"RecoCosineZ",kRecoCosZ},
    {"OscillationChannel",kOscChannel},
    {"Mode",kMode}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kTrueCosZ,"TrueCosineZ"},    
    {kRecoCosZ,"RecoCosineZ"},
    {kOscChannel,"OscillationChannel"},
    {kMode,"Mode"}
  };
  
  /// Array filled with MC samples for each oscillation channel
  std::vector<struct dunemc_base> dunemcSamples;

  /// Is the sample e-like
  bool IsELike;

  /// Multiplicative scaling to scale from the assumed 400ktyr value in the CAF files
  double ExposureScaling;
};

#endif
