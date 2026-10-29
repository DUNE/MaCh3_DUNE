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
    kRecoNeutrinoEnergy,
    kTrueCosZ,
    kRecoCosZ,
    kOscChannel,
    kMode,
    kTargetNucleus,
    kMinDistToWall,
    kRecoInelasticity
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
    {"TrueCosineZ",kTrueCosZ},
    {"RecoCosineZ",kRecoCosZ},
    {"OscillationChannel",kOscChannel},
    {"Mode",kMode},
    {"TargetNucleus", kTargetNucleus},
    {"MinDistToWall", kMinDistToWall},
    {"RecoInelasticity", kRecoInelasticity}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kRecoHadronEnergy,"RecoHadronEnergy"},
    {kRecoLeptonEnergy,"RecoLeptonEnergy"},
    {kTrueCosZ,"TrueCosineZ"},    
    {kRecoCosZ,"RecoCosineZ"},
    {kOscChannel,"OscillationChannel"},
    {kMode,"Mode"},
    {kTargetNucleus, "TargetNucleus"},
    {kMinDistToWall, "MinDistToWall"},
    {kRecoInelasticity, "RecoInelasticity"}
  };
  
  /// Array filled with MC samples for each oscillation channel
  std::vector<dunemc_atm> dunemcSamples;

  /// Is the sample e-like
  std::vector<int> IsELike;

  /// Multiplicative scaling to scale from the assumed 400ktyr value in the CAF files
  double ExposureScaling;

  /// Path to the input spline file (for per-event spline mode)
  std::vector<std::string> fInputSplines;
  
  /// Build the per-event splines from the xsec variations stored in the CAFs (SRTrueInteraction::xsec_systs)
  bool fSplinesFromCAF;

  /// Knot layout of one CAF xsec dial used as a spline parameter
  struct CAFSplineDial {
    /// Position of the dial in SRWeightGlobal::xsec_params and SRTrueInteraction::xsec_systs
    size_t CAFIndex;
    /// Number of variations stored in the CAF for this dial
    size_t nCAFVariations;
    /// For each sorted knot, the CAF variation it comes from (-1 for an added nominal knot with weight 1)
    std::vector<int> KnotSource;
    /// Sorted knot positions
    std::vector<double> Knots;
    /// Weight capping from the spline parameter config
    double KnotLowBound;
    double KnotUpBound;
  };
  std::vector<CAFSplineDial> fCAFSplineDials;
  /// Names and variations of the CAF xsec dials, used to check every input file has the same dials
  std::vector<std::pair<std::string, std::vector<float>>> fCAFDialHeaders;
  /// Per-event splines built while reading the CAFs, handed over to MonolithSplineHandlerDUNE in SetupSplines
  std::vector<std::vector<TResponseFunction_red*>> fCAFSplines;
  /// Memory of the CAF splines (objects and knots), released once MonolithSplineHandlerDUNE has packed them
  SplineArena fCAFSplineArena;

  /// Read the xsec dial headers from the globalTree of a CAF file
  std::vector<std::pair<std::string, std::vector<float>>> ReadCAFDialHeaders(TFile* File) const;
  /// Match the spline parameters of this sample to the CAF xsec dials
  void SetupCAFSplineDials(TFile* File);
  
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
  
  int ReturnSampleIdentifier(const std::vector<double>& CVNScores, double MinDistanceToWall);
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
  enum Variables{neutrinoenergy_reco,hadronicenergy_reco,leptonicenergy_reco,coszenith_reco,SampleIndex,neutrinoflavour_interacting,neutrinoflavour_generated,OscChannelIndex,mode,IsCC,Target,neutrinoenergy_true,coszenith_true,weights_flux,MinDistToWall,eventID,nVariables};

  std::string InputFileName;
  std::string InputFileDirectory;  
  
  /// File path to binary file which contains Eigen::MatrixXd
  std::string EigenInputFile;

  std::vector<std::vector<std::string>> BinningVariables;
};

#endif
