#ifndef _SampleHandlerAtm_h_
#define _SampleHandlerAtm_h_

#include "Splines/SplineHandlerFactoryDUNE.h"
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
    kRecoHadronEnergy,
    kRecoLeptonEnergy,
    kRecoY,
    kRecoNeutrinoEnergy,
    kTrueCosZ,
    kRecoCosZ,
    kOscChannel,
    kMode,
    kTargetNucleus,
    kMinDistToWall
  };

  /// @brief Convert dunemc_atm struct object to Eigen::MatrixXd and write to binary file
  void TransferToEigen(std::string FileName);

  BinningHandler* GetBinningHandler() const {
    return Binning.get();
  }

  std::vector<std::string> GetSampleBinningVariables(int iSample) {
    return BinningVariables.at(iSample);
  }

protected:
  /// @brief Initialises object
  void Init();

  /// @brief Initialise data hist (can be overridden)
  void InititialiseData() override;

  /// @brief Function to setup MC from file
  /// @param iSample sample ID
  /// @return Total number of events
  int SetupExperimentMC();

  /// @brief Read binary file into Eigen::MatrixXd and transfer into dunemc_atm struct
  int ReadFromEigen();

  /// @brief Tells FD base which variables to point to/be set to
  /// @param iSample Sample ID
  void SetupMC();

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void AddAdditionalWeightPointers();

  /// @brief Sets up splines 
  void SetupSplines();

  /// @brief Initialise per-event spline weight pointers (for MonolithSplineHandler)
  void InitialiseSplineObjectPerEvent();

  /// @brief Cleanup memory
  void CleanMemoryBeforeFit() override {};

  //=============================================================================
  // Fast atmospheric reweight (FastOptions)
  //=============================================================================
public:
  /// @brief Fast reweight: same weights and histograms as SampleHandlerBase::Reweight, computed from compact
  ///        per-event arrays built once (selection and bin are fixed because no kinematic shifts exist here)
  void Reweight() override;
protected:
  /// @brief Keeps GetEventWeight (plots, projections) correct when the fast reweight is used: recomputes this
  ///        event's detector and flux systematic weights from the current parameter values
  void CalcWeightFunc(const int iEvent) override;
  /// @brief Read the fast-build options (FastOptions node, InputFiles:SplineCacheFile, M3_* environment overrides)
  void ReadFastOptions();
  /// @brief Build the compact arrays used by the fast reweight
  void BuildFastEventArrays();
  /// @brief Fill SampleHandler_array from the compact arrays
  void FillArrayFast();
  /// @brief Whether the functional shift with these parameters applies to event iEvent (same test as core registration)
  bool FunctionalShiftApplies(const std::vector<std::string>& ParNames, const int iEvent);

  bool UseFastReweight = true;
  bool UseFastSplines = true;
  bool FastArraysReady = false;
  FastAtmSplineOptions FastSplineOptions;
  std::string InputSplinesMD5Sum;
  std::vector<int> DetectorSystematicParameterIndices;
  std::vector<int> AtmFluxParameterIndices;
  /// Per event: bit 0 = detector shift applies, bit 1 = flux shift applies (same test as core registration)
  std::vector<uint8_t> FunctionalApplies;

  /// Compact per-event information (only events that can enter a histogram)
  struct FastEventArrays {
    std::vector<int> Bin;                    ///< global bin
    std::vector<double> Const;               ///< product of the constant weight pointers (flux weight, exposure)
    std::vector<uint32_t> DynOff;            ///< offsets into DynPtr (size N+1)
    std::vector<const M3::float_t*> DynPtr;  ///< changing weights: oscillation, splines, norms
    std::vector<double> DetRm1;              ///< [iSyst*N + e]: detector ratio - 1 (0 if the shift does not apply)
    std::vector<double> FluxFm1;             ///< flux factor - 1 (0 if none)
    std::vector<uint8_t> FluxPar;            ///< index into AtmFluxParameterIndices, 255 = none
    size_t N = 0;
    size_t Bytes() const;
  } Fast;


  void RegisterFunctionalParameters();
  
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
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(const int KinematicVariable, const int iEvent) const override;

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const override;

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
    {"RecoHadronEnergy",kRecoHadronEnergy},
    {"RecoLeptonEnergy",kRecoLeptonEnergy},
    {"RecoY",kRecoY},    
    {"TrueCosineZ",kTrueCosZ},
    {"RecoCosineZ",kRecoCosZ},
    {"OscillationChannel",kOscChannel},
    {"Mode",kMode},
    {"TargetNucleus", kTargetNucleus},
    {"MinDistToWall", kMinDistToWall}    
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kRecoHadronEnergy,"RecoHadronEnergy"},
    {kRecoLeptonEnergy,"RecoLeptonEnergy"},
    {kRecoY,"RecoY"},
    {kTrueCosZ,"TrueCosineZ"},    
    {kRecoCosZ,"RecoCosineZ"},
    {kOscChannel,"OscillationChannel"},
    {kMode,"Mode"},
    {kTargetNucleus, "TargetNucleus"},
    {kMinDistToWall, "MinDistToWall"}
  };
  
  /// Array filled with MC samples for each oscillation channel
  std::vector<dunemc_atm> dunemcSamples;

  /// Is the sample e-like
  std::vector<int> IsELike;

  /// Multiplicative scaling to scale from the assumed 400ktyr value in the CAF files
  double ExposureScaling;

  /// Path to the input spline file (for per-event spline mode)
  std::vector<std::string> fInputSplines;
  
  /// Enums to define event selections
  enum EventSelectionIndices {
    kEventSel_Unknown = -1,    
    kEventSel_FC_NuE,
    kEventSel_FC_NuMu,
    kEventSel_FC_NC,
    kEventSel_PC_NuE,
    kEventSel_PC_NuMu,
    kEventSel_PC_NC,    
    nEventSelections,
  };

  /// Enums to define ordering of the CVN scores from the CAF files
  enum CVNScoreIndices {
    kCVN_NuE,
    kCVN_NuMu,
    kCVN_NC,
    nCVN_Scores
  };

  /// Cut value to separate Fully Contained and Partially Contained events
  double FCPCSeparation;
  
  int ReturnSampleIdentifier(std::vector<double> CVNScores, double MinDistanceToWall);
  std::vector<std::string> EventSelectionNames = std::vector<std::string>(nEventSelections);
  std::vector<int> EventSelection_to_SampleIndex_Map = std::vector<int>(nEventSelections,kEventSel_Unknown);

  void SetupDetectorSystematicRatios();
  std::string fDetectorSystematicsFileName;
  std::vector<std::string> DetectorSystematicParameterNames;
  std::unordered_map<std::string,std::vector<int>> DetectorSystematicDimensionMap;

  void SetupAtmFluxSplineSystematics();
  std::string fHONDAFluxUncerSplineFileName;
  double fHONDAFluxUncerSeparationValue;
  std::vector<std::string> AtmFluxParameterNames;  
  
  /// Enum defining CAF->Eigen transformation
  enum Variables{neutrinoenergy_reco,hadronicenergy_reco,leptonicenergy_reco,inelasticity_reco,coszenith_reco,SampleIndex,neutrinoflavour_interacting,neutrinoflavour_generated,OscChannelIndex,mode,IsCC,Target,neutrinoenergy_true,coszenith_true,weights_flux,MinDistToWall,eventID,nVariables};

  std::string InputFileName;
  std::string InputFileDirectory;  
  
  /// File path to binary file which contains Eigen::MatrixXd
  std::string EigenInputFile;

  std::vector<std::vector<std::string>> BinningVariables;
};

#endif
