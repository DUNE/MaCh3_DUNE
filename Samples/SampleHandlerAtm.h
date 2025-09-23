#ifndef _SampleHandlerAtm_h_
#define _SampleHandlerAtm_h_

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerFD.h"

#include "StructsDUNE.h"
/// @brief Base class for handling atmospheric samples
class SampleHandlerAtm : virtual public SampleHandlerFD // Inherits from SampleHandlerFD
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
  // Define available kinetmatic variables that can be used for binning and cuts
  enum KinematicTypes{kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueCosZ,kRecoCosZ,kOscChannel,kMode};
  
protected: // Core functions
  /// @brief Initialises object
  void Init();

  /// @brief Function to setup MC from file
  /// @param iSample sample ID
  /// @return Total number of events
  int SetupExperimentMC();

  /// @brief Tells FD base which variables to point to/be set to
  /// @param iSample Sample ID
  void SetupFDMC();

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  // void SetupWeightPointers();
  void SetupWeightPointers() override; // muyuan changed this. without voerride, the compiler
  //treats it as a new function rather than overriding the base class function, so when the 
  /// base class function is called, it doesn't call this one

  /// @brief Sets up splines 
  void SetupSplines();

  /// @brief Cleanup memory
  void CleanMemoryBeforeFit() override {};


  /// @brief Muyuan's function to register functional parameters
  /// @details Not implemented, but required by the base class
  // void RegisterFunctionalParameters() override {};
  void RegisterFunctionalParameters() override;
  
  //DB functions which could be initialised to do something which is non-trivial
  
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iSample sample ID
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iEvent) {(void) iEvent; return 1.;} // always return 1.0 (no effect)
  
  /// @brief NOT IMPLEMENTED: Apply kinematic shifts
  /// @param iSample Sample Number
  /// @param iEvent Event number
  void applyShifts(int iEvent) {(void)iEvent;}
  
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinPar Kinematic parameter enum val
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(int KinematicVariable, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent);

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
  std::vector<dunemc_base> dunemcSamples;

  /// Is the sample e-like
  bool IsELike;

  /// Multiplicative scaling to scale from the assumed 400ktyr value in the CAF files
  double ExposureScaling;


private: // MUYUAN made

  struct SplineSet {
    TSpline3* energy_spline;
    TSpline3* cosZ_spline;
    double energy_min, energy_max;
    double cosZ_min, cosZ_max;
  };
  SplineSet spline_sets[4]; // [0]=nue, [1]=anue, [2]=numu, [3]=anumu
  // MUYUAN: Test spline for atmospheric flux
  // TSpline3* spline_nue_E;
  // TSpline3* spline_nue_Cos;
  // TSpline3* spline_anue_E;
  // TSpline3* spline_anue_Cos;
  // TSpline3* spline_numu_E;
  // TSpline3* spline_numu_Cos;
  // TSpline3* spline_anumu_E;
  // TSpline3* spline_anumu_Cos;

  // MUYUAN: added to SetupExperimentMC, info from StructsDUNE.h
  std::vector<double> original_flux_weights;

  void AtmFluxShift(const double* par, std::size_t iSample, std::size_t iEvent);


  // MUYUAN: Parameter values for flux variations
  double param_atmflux_nue;
  double param_atmflux_anue;
  double param_atmflux_numu;
  double param_atmflux_anumu;
  
  // MUYUAN: Spline loading functions
  void LoadFluxSplines();
  TSpline3* LoadSingleSpline(const std::string& filename, const std::string& splinename);
  
  // MUYUAN: Flux weight calculation
  double CalculateFluxWeight(int iEvent);

};











#endif
