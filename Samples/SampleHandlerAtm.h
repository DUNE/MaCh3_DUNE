#ifndef _SampleHandlerAtm_h_
#define _SampleHandlerAtm_h_

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerFD.h"
#include "StructsDUNE.h"

#include <array>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

// Forward declarations
class TSpline3;
class TH2D;

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
  // Define available kinematic variables that can be used for binning and cuts
  enum KinematicTypes{kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueCosZ,kRecoCosZ,kOscChannel,kMode};

protected: // Core functions
  /// @brief Initialises object
  void Init();

  /// @brief Function to setup MC from file
  /// @return Total number of events
  int SetupExperimentMC();

  /// @brief Tells FD base which variables to point to/be set to
  void SetupFDMC();

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void SetupWeightPointers() override;

  /// @brief Sets up splines 
  void SetupSplines();

  /// @brief Cleanup memory
  void CleanMemoryBeforeFit() override {};

  /// @brief Register functional parameters
  void RegisterFunctionalParameters() override;
  
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iEvent) {(void) iEvent; return 1.;} // always return 1.0 (no effect)
  
  /// @brief NOT IMPLEMENTED: Apply kinematic shifts
  /// @param iEvent Event number
  void applyShifts(int iEvent) {(void)iEvent;}
  
  /// @brief Returns pointer to kinematic parameter for event in Structs DUNE
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) override;
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) override;

  /// @brief Returns kinematic parameter value for event
  double ReturnKinematicParameter(int KinematicVariable, int iEvent) override;
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent) override;

  /// @brief Gets binning for a given parameter
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameterStr);
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

private:
  // OPTIMIZED: Spline storage structure
  struct SplineSet {
    TSpline3* energy_spline;
    TSpline3* cosZ_spline;
    double energy_min, energy_max;
    double cosZ_min, cosZ_max;
  };
  
  // OPTIMIZED: Pre-computed spline lookup grids (much faster than TSpline3::Eval)
  struct SplineGrid {
    TH2D* energy_cosZ_grid;  // 2D histogram for fast lookup
    double energy_min, energy_max;
    double cosZ_min, cosZ_max;
    int energy_bins, cosZ_bins;
    double energy_step, cosZ_step;
  };
  
  SplineSet spline_sets[4];   // [0]=nue, [1]=anue, [2]=numu, [3]=anumu
  SplineGrid spline_grids[4]; // [0]=nue, [1]=anue, [2]=numu, [3]=anumu
  
  // Cache for event spline responses (computed once, reused)
  std::vector<std::array<double, 4>> cached_spline_responses; // [event][neutrino_type]
  bool responses_cached;

  // Store original flux weights for restoration
  std::vector<double> original_flux_weights;

  // Flux parameters
  double param_atmflux_nue;
  double param_atmflux_anue;
  double param_atmflux_numu;
  double param_atmflux_anumu;
  
  // Helper methods
  void LoadFluxSplines();
  void PrecomputeSplineGrids();
  void CacheAllSplineResponses();
  double GetCachedFluxWeight(int iEvent, int pdg_index) const;
  double CalculateFluxWeight(int iEvent);
  void AtmFluxShift(const double* par, std::size_t iSample, std::size_t iEvent);
};

#endif