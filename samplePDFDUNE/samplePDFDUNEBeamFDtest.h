#ifndef _samplePDFDUNEBeamFD_h_
#define _samplePDFDUNEBeamFD_h_

#include "splines/splinesDUNE.h"
#include "samplePDF/samplePDFFDBase.h"
#include <vector>
#include <string>


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

  std::unordered_map<std::string, std::pair<int, FuncParFuncType>> functionalParameterRegistry;


  /// @brief Enum to identify kinematics
  enum KinematicTypes {kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueXPos,kTrueYPos,kTrueZPos,kCVNNumu,kCVNNue,kM3Mode,kOscChannel,kIsFHC, kq0, kq3, k_pT, k_pz, k_global_bin_number, kp_lep, ktheta_lep, kELepRec, kEHadRec, kERec_minus_Etrue, kERecQE, kENuProxy_minus_Enutrue, kyRec, keHad_av, kisCC, kisRelativeEnubias};
  std::pair<std::vector<double>, std::vector<double>> Return2DKinematicParameterBinning(std::string KinematicParameterStr);
  std::vector<double> f1DEdges;
  TH1D* f1DHist = nullptr;
  void Setup1DHist(); // optional helper to build fQ0Q3Hist
  TH1D* Get1DHist() const { return f1DHist; }


  // Global energy scale systematics
  void TotalEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  
  //void TotalEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  
  //void TotalEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  
  // Particle specific energy uncertainties
  // Charged hadron
  //void HadEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  //void HadEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  //void HadEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  // Muons
  //void MuEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  //void MuEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  //void MuEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  // Neutrons
  //void NEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  //void NEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  //void NEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  // Electromagnetic showers
  //void EMEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  
  //void EMEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  
  //void EMEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent);
  

  // Resolution uncertainties
  //void HadRes(const double * par, std::size_t iSample, std::size_t iEvent);
  //void MuRes(const double * par, std::size_t iSample, std::size_t iEvent);
  //void NRes(const double * par, std::size_t iSample, std::size_t iEvent);
  //void EMRes(const double * par, std::size_t iSample, std::size_t iEvent);
  //void EMResCCNue(const double * par, std::size_t iSample, std::size_t iEvent);

  
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

  void RegisterFunctionalParameters();
  void RegisterIndividualFunctionalParameter(const std::string& fpName, int fpEnum, FuncParFuncType fpFunc);
  double CalculatePOT();
  
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
 
  // dunemc
  /// DUNE MC sampels
  std::vector<struct dunemc_base> dunemcSamples;

  /// Value of POT used for sample
  double pot;
  bool iselike;
  double isFHC;

  //Positions of FD Detector systematics
  double tot_escale_nd_pos;
  double tot_escale_sqrt_nd_pos;
  double tot_escale_invsqrt_nd_pos;
  double had_escale_nd_pos;
  double had_escale_sqrt_nd_pos;
  double had_escale_invsqrt_nd_pos;
  double mu_escale_nd_pos;
  double mu_escale_sqrt_nd_pos;
  double mu_escale_invsqrt_nd_pos;
  double n_escale_nd_pos;
  double n_escale_sqrt_nd_pos;
  double n_escale_invsqrt_nd_pos;
  double em_escale_nd_pos;
  double em_escale_sqrt_nd_pos;
  double em_escale_invsqrt_nd_pos;
  double had_res_nd_pos;
  double mu_res_nd_pos;
  double n_res_nd_pos;
  double em_res_nd_pos;
  double cvn_numu_nd_pos;
  double cvn_nue_nd_pos;

  /// FD Detector Systematics
  std::vector<const double*> FDDetectorSystPointers;

  /// Number of FD Detector Systematics
  int nFDDetectorSystPointers;

  // === HH: Functional parameters ===
  enum FuncParEnum {kTotalEScale,
    kTotalEScaleSqrt,
    kTotalEScaleInvSqrt, 
    kHadEScale, kHadEScaleSqrt, kHadEScaleInvSqrt,
    kMuEScale, kMuEScaleSqrt, kMuEScaleInvSqrt,
    kNEScale, kNEScaleSqrt, kNEScaleInvSqrt,
    kEMEScale, 
    kEMEScaleSqrt, 
    kEMEScaleInvSqrt,
    kHadRes, kMuRes, kNRes, kEMRes
  };

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
    {"q0",kq0},
    {"q3",kq3},
    {"pT",k_pT},
    {"pz",k_pz},
    {"global_bin_number",k_global_bin_number},
    {"p_lep",kp_lep},
    {"theta_lep",ktheta_lep},
    {"ELepRec",kELepRec},
    {"EHadRec",kEHadRec},
    {"ERec_minus_Etrue",kERec_minus_Etrue},
    {"ERecQE",kERecQE},
    {"ENuProxy_minus_Enutrue",kENuProxy_minus_Enutrue},
    {"yRec",kyRec},
    {"eHad_av",keHad_av},
    {"isCC", kisCC},
    {"isRelativeEnubias",  kisRelativeEnubias}
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
    {kq0, "q0"},
    {kq3, "q3"},
    {k_pT, "pT"},
    {k_pz, "pz"},
    {k_global_bin_number, "global_bin_number"},
    {kp_lep, "p_lep"},
    {ktheta_lep, "theta_lep"},
    {kELepRec, "ELepRec"},
    {kEHadRec,"EHadRec"},
    {kERec_minus_Etrue, "ERec_minus_Etrue"},
    {kERecQE, "ERecQE"},
    {kENuProxy_minus_Enutrue, "ENuProxy_minus_Enutrue"},
    {kyRec, "yRec"},
    {keHad_av, "eHad_av"},
    {kisCC, "isCC"},
    {kisRelativeEnubias, "isRelativeEnubias"}
  };
};



#endif
