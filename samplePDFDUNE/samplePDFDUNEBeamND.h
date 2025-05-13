#ifndef _samplePDFDUNEBeamND_h_
#define _samplePDFDUNEBeamND_h_

#include "splines/splinesDUNE.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"
#include "samplePDF/Structs.h"

class samplePDFDUNEBeamND : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamND(std::string mc_version, covarianceXsec* xsec_cov, TMatrixD* nd_cov, covarianceOsc* osc_cov) ;
  ~samplePDFDUNEBeamND();

  enum KinematicTypes {kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kyRec,kOscChannel,kMode,kIsFHC,
    kIsCC,kRecoNumu,kRecoNue,kCCNumu,kCCNue,kNotCCNumu,kRecoQ,
    kRecoHadEnergy,kRecoLepEnergy,kRecoPEnergy,
    kRecoPipEnergy,kRecoPimEnergy,kRecoPi0Energy,kRecoNEnergy, kNuPDG,
  };
  
 protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();

  // === HH: Functional parameters ===
  enum FuncParEnum {kDebugNothing, kDebugShift, 
    kTotalEScale, kTotalEScaleNotCCNumu, 
    kTotalEScaleSqrt, kTotalEScaleSqrtNotCCNumu, 
    kTotalEScaleInvSqrt, kTotalEScaleInvSqrtNotCCNumu,
    kHadEScale, kHadEScaleSqrt, kHadEScaleInvSqrt,
    kMuEScale, kMuEScaleSqrt, kMuEScaleInvSqrt,
    kNEScale, kNEScaleSqrt, kNEScaleInvSqrt,
    kEMEScale, kEMEScaleCCNue, 
    kEMEScaleSqrt, kEMEScaleSqrtCCNue,
    kEMEScaleInvSqrt, kEMEScaleInvSqrtCCNue,
    kHadRes, kMuRes, kNRes, kEMRes, kEMResCCNue
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

  // Debugging
  void DebugShift(const double * par, std::size_t iSample, std::size_t iEvent);
  // =================================
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  double ReturnKinematicParameter(int KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
  double ReturnKinematicParameter(KinematicTypes Kinpar, int iSample, int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.; (void)iSample; (void)iEvent;}

  // HH - Replacing poisson likelihood with gaussian for ND
  double GetLikelihood() override;
  TMatrixD *NDCovMatrix;
  std::vector<std::vector<double>> NDInvertCovMatrix;
  std::vector<double> FlatDataMCDiff;
  void setNDCovMatrix();
  bool isNDCovSet = false;
  int nXBins, nYBins, covSize;

  std::vector<struct dunemc_base> dunendmcSamples;

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
    {"yRec",kyRec},
    {"OscillationChannel",kOscChannel},
    {"Mode",kMode},
    {"IsFHC",kIsFHC}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
    {kyRec,"yRec"},
    {kOscChannel,"OscillationChannel"},
    {kMode,"Mode"},
    {kIsFHC,"IsFHC"}
  };

  TString _nutype;
  int _mode;

  double pot;

  // dunendmc Variables
  double _ev;
  double _erec;
  double _erec_lep;
  double _erec_had;
  int _reco_numu;
  int _reco_nue;

  double _eRecoP;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoN;

  double _LepNuAngle;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _BeRPA_cvwgt;
  int _isCC;
  int _nuPDGunosc;
  int _nuPDG;
  int _run;
  int _isND;
  int _isFHC;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;
  double _LepTheta;
  double _Q2;
  int _reco_q;

  // configuration 
  bool iselike;
  bool isND;
  double IsFHC;

};

#endif
