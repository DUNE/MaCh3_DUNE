#ifndef _samplePDFDUNEBeamFD_h_
#define _samplePDFDUNEBeamFD_h_

#include "splines/splinesDUNE.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"
/// @brief Base class for handling FD Beam samples
class samplePDFDUNEBeamFD : virtual public samplePDFFDBase
{
public:


  /// @brief samplePDFDUNE FD beam Constructor
  /// @param mc_version Config Name
  /// @param xsec_cov Cross-section covariance matrix
  /// @param osc_cov Oscillation covariance matrix
  samplePDFDUNEBeamFD(std::string mc_version, covarianceXsec* xsec_cov, covarianceOsc* osc_cov);

  /// @brief destructor
  ~samplePDFDUNEBeamFD();

  /// @brief Enum to identify kinematics
  enum KinematicTypes {kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueXPos,kTrueYPos,kTrueZPos,kCVNNumu,kCVNNue,kM3Mode,kOscChannel,kIsFHC, kTrueCCnue, kTrueCCnumu};

protected:
  /// @brief Initialises object
  void Init();

  /// @brief Function to setup MC from file
  /// @param iSample sample ID
  /// @return Total number of events
  int setupExperimentMC(int iSample);

  /// @brief Tells FD base which variables to point to/be set to
  /// @param iSample Sample ID
  void setupFDMC(int iSample);

  /// @brief Sets up pointers weights for each event (oscillation/xsec/etc.)
  void SetupWeightPointers();
  void SetupSplines();
	
	// === HH: Functional parameters ===
  enum FuncParEnum {kTotalEScale, kTotalEScaleNotCCNumu, 
    kTotalEScaleSqrt, kTotalEScaleSqrtNotCCNumu, 
    kTotalEScaleInvSqrt, kTotalEScaleInvSqrtNotCCNumu,
    kHadEScale, kHadEScaleSqrt, kHadEScaleInvSqrt,
    kMuEScale, kMuEScaleSqrt, kMuEScaleInvSqrt,
    kNEScale, kNEScaleSqrt, kNEScaleInvSqrt,
    kEMEScale, kEMEScaleCCNue, 
    kEMEScaleSqrt, kEMEScaleSqrtCCNue,
    kEMEScaleInvSqrt, kEMEScaleInvSqrtCCNue,
    kHadRes, kMuRes, kNRes, kEMRes, kEMResCCNue,
	kRecoCVNNumu, kRecoCVNNue
  };
  void RegisterFunctionalParameters();
  void resetShifts(int iSample, int iEvent);

  // Global energy scale systematics
  void TotalEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  void TotalEScaleNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent);
  void TotalEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  void TotalEScaleSqrtNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent);
  void TotalEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  void TotalEScaleInvSqrtNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent);

  // Particle specific energy uncertainties
  // Charged hadron
  void HadEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  void HadEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  void HadEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  // Muons
  void MuEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  void MuEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  void MuEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  // Neutrons
  void NEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  void NEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  void NEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  // Electromagnetic showers
  void EMEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  void EMEScaleCCNue(const double * par, std::size_t iSample, std::size_t iEvent);
  void EMEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  void EMEScaleSqrtCCNue(const double * par, std::size_t iSample, std::size_t iEvent);
  void EMEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  void EMEScaleInvSqrtCCNue(const double * par, std::size_t iSample, std::size_t iEvent);

  // Resolution uncertainties
  void HadRes(const double * par, std::size_t iSample, std::size_t iEvent);
  void MuRes(const double * par, std::size_t iSample, std::size_t iEvent);
  void NRes(const double * par, std::size_t iSample, std::size_t iEvent);
  void EMRes(const double * par, std::size_t iSample, std::size_t iEvent);
  void EMResCCNue(const double * par, std::size_t iSample, std::size_t iEvent);

  //Reconstruction (CVN) uncertainties
  void RecoCVNNumu(const double * par, std::size_t iSample, std::size_t iEvent);
  void RecoCVNNue(const double * par, std::size_t iSample, std::size_t iEvent);

  /*
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter Type
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event 
  double ReturnKinematicParameter (KinematicTypes KinPar, int iSample, int iEvent);*/

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as int
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event 
  double ReturnKinematicParameter (int KinematicVariable, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinPar Kinematic Parameter Type
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent); 

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  inline std::string ReturnStringFromKinematicParameter(int KinematicParameterStr);
  
  //DB functions which could be initialised to do something which is non-trivial

  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iSample sample ID
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iSample, int iEvent) {(void) iSample; (void)iEvent; return 1.;}

  /// @brief Apply kinematic shifts
  /// @param iSample Sample Number
  /// @param iEvent Event number
  void applyShifts(int iSample, int iEvent);

  // dunemc
  /// DUNE MC sampels
  std::vector<struct dunemc_base> dunemcSamples;

  /// Value of POT used for sample
  double pot;
  bool iselike;
  double isFHC;

  //Positions of FD Detector systematics
  double tot_escale_fd_pos;
  double tot_escale_sqrt_fd_pos;
  double tot_escale_invsqrt_fd_pos;
  double had_escale_fd_pos;
  double had_escale_sqrt_fd_pos;
  double had_escale_invsqrt_fd_pos;
  double mu_escale_fd_pos;
  double mu_escale_sqrt_fd_pos;
  double mu_escale_invsqrt_fd_pos;
  double n_escale_fd_pos;
  double n_escale_sqrt_fd_pos;
  double n_escale_invsqrt_fd_pos;
  double em_escale_fd_pos;
  double em_escale_sqrt_fd_pos;
  double em_escale_invsqrt_fd_pos;
  double had_res_fd_pos;
  double mu_res_fd_pos;
  double n_res_fd_pos;
  double em_res_fd_pos;
  double cvn_numu_fd_pos;
  double cvn_nue_fd_pos;

  /// FD Detector Systematics
  std::vector<const double*> FDDetectorSystPointers;

  /// Number of FD Detector Systematics
  int nFDDetectorSystPointers;

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
	{kTrueCCnue,"IsTrueCCnue"},
	{kTrueCCnumu,"IsTrueCCnumu"}
  };
};

#endif
