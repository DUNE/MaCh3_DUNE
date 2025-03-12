#ifndef _samplePDFDUNEBeamND_h_
#define _samplePDFDUNEBeamND_h_

#include "splines/splinesDUNE.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"

class samplePDFDUNEBeamND : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamND(std::string mc_version, covarianceXsec* xsec_cov, TMatrixD* nd_cov, covarianceOsc* osc_cov) ;
  ~samplePDFDUNEBeamND();

  enum KinematicTypes {kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kyRec,kOscChannel,kMode,kIsFHC};
  
 protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.; (void)iSample; (void)iEvent;}
  void applyShifts(int iSample, int iEvent);

  void setNDCovMatrix();
  double GetLikelihood();

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


  //TFile *_sampleFile;
  //TTree *_data;

  TChain* _data;
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

  //Positions of ND Detector systematics
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
  // The ND detector covariance matrix
  //
  TMatrixD *NDCovMatrix;
  
  // The inverse ND detector covariance matrix
  TMatrixD *NDInvCovMatrix;
  bool isNDCovSet = false;

  double **NDInvertCovMatrix;

  std::vector<const double*> NDDetectorSystPointers;
  int nNDDetectorSystPointers;
};

#endif
