#pragma once

#include "Splines/BinnedSplineHandlerDUNE.h"

#include "Samples/SampleHandlerFD.h"

#include "Samples/StructsDUNE.h"
/// @brief Base class for handling FD Beam samples
class SampleHandlerBeamOffAxis : virtual public SampleHandlerFD {
public:
  /// @brief SampleHandler FD beam Constructor
  /// @param mc_version Config Name
  /// @param xsec_cov Cross-section covariance matrix
  /// @param osc_cov Oscillation covariance matrix
  /// @param Oscillator_ Shared Oscillation Handler object
  SampleHandlerBeamOffAxis(
      std::string mc_version, ParameterHandlerGeneric *xsec_cov,
      const std::shared_ptr<OscillationHandler> &Oscillator_);

  /// @brief destructor
  ~SampleHandlerBeamOffAxis();

  // below is ugly, but lets us define it only once and get the enum and both
  // maps https://en.wikipedia.org/wiki/X_macro

#define LIST_OF_VARIABLES                                                      \
  X(TrueNeutrinoEnergy)                                                        \
  X(RecoNeutrinoEnergy)                                                        \
  X(TrueXPos)                                                                  \
  X(TrueYPos)                                                                  \
  X(TrueZPos)                                                                  \
  X(M3Mode)                                                                      \
  X(IsFHC)                                                                     \
  X(ELepRec)                                                                   \
  X(Enubias)

#define X(a) k##a,

  /// @brief Enum to identify kinematics
  enum KinematicTypes { LIST_OF_VARIABLES };

#undef X

protected:
  /// @brief Initialises object
  void Init();

  /// @brief Function to setup MC from file
  /// @return Total number of events
  int SetupExperimentMC();

  /// @brief Tells FD base which variables to point to/be set to
  void SetupFDMC();

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void SetupWeightPointers();
  void SetupSplines();

  enum FuncParEnum {
    kTotalEScaleND = 0,
    kTotalEScaleND_invsqrt,
    kTotalEScaleND_sqrt,
    kTotalEScaleND_mu,
    kTotalEScaleND_musqrt,
    kTotalEScaleND_muinvsqrt,
    kTotalEScaleND_had,
    kTotalEScaleND_hadsqrt,
    kTotalEScaleND_hadinvsqrt,
    kTotalEScaleND_EM,
    kTotalEScaleND_EMinvsqrt,
    kTotalEScaleND_EMsqrt,
    kMuonRes_ND,
    kNRes_ND,
    kHadRes_ND,
    kNFuncPars
  };

  void RegisterFunctionalParameters() override;
  void resetShifts(int iEvent) override;

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter Type
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(KinematicTypes KinPar, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as int
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(int KinematicVariable, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast ->
  /// int)
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinPar Kinematic Parameter Type
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double *GetPointerToKinematicParameter(KinematicTypes KinPar,
                                               int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast ->
  /// int)
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double *GetPointerToKinematicParameter(std::string KinematicParameter,
                                               int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter as double (gets cast -> int)
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double *GetPointerToKinematicParameter(double KinematicVariable,
                                               int iEvent);

  // std::vector<double> ReturnKinematicParameterBinning(std::string
  // KinematicParameter);
  inline std::string
  ReturnStringFromKinematicParameter(int KinematicParameterStr);

  // DB functions which could be initialised to do something which is
  // non-trivial

  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iEvent) {
    (void)iEvent;
    return 1.;
  }

  std::vector<std::vector<std::vector<std::vector<TH2D*>>>> GetBinnedWeights(std::vector<std::string> ParamNames, std::vector<std::vector<int>> ParamModes, std::vector<double> TrueEBins);  

  std::vector<dunemc_beamoffaxis> dunemcSamples;

  /// Value of POT used for sample
  double pot;
  double isFHC;

#define X(a) {#a, k##a},
  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
      LIST_OF_VARIABLES};

#undef X
#define X(a) {k##a, #a},
  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
      LIST_OF_VARIABLES};

#undef X
#undef LIST_OF_VARIABLES

  /// @brief Cleanup memory
  void CleanMemoryBeforeFit() override {};
};
