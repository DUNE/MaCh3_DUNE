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

  enum KinematicTypes {
    kTrueNeutrinoEnergy, 
    kRecoQ, 
    kRecoNeutrinoEnergy, 
    kIsFHC,
    kyRec, 
    kOscChannel, 
    kMode, 
    kMuonMom, 
    kMuonEnergy, 
    kRecoMuonEnergy, 
    kMuonTheta, 
    kPipMom,
    kPipEnergy,
    kRecoPipEnergy, 
    kPipTheta, 
    kMuonEDiff, 
    kPipEDiff, 
    kEDiff,
    kStartX,
    kStartY,
    kStartZ
  };  
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
  double GetLikelihood() override;

  std::vector<struct dunemc_base> dunendmcSamples;

  const std::unordered_map<std::string, KinematicTypes> KinematicParameterMap = {
    {"yRec", kyRec},
    {"TrueNeutrinoEnergy", kTrueNeutrinoEnergy},
    {"RecoQ", kRecoQ},
    {"RecoNeutrinoEnergy", kRecoNeutrinoEnergy},
    {"IsFHC", kIsFHC},
    {"OscChannel", kOscChannel},
    {"Mode", kMode},
    {"MuonMom", kMuonMom},
    {"MuonEnergy", kMuonEnergy},
    {"MuonTheta", kMuonTheta},
    {"RecoMuonEnergy", kRecoMuonEnergy},
    {"PipMom", kPipMom},
    {"PipEnergy", kPipEnergy},
    {"RecoPipEnergy", kRecoPipEnergy},
    {"PipTheta", kPipTheta},
    {"MuonEDiff", kMuonEDiff},
    {"PipEDiff", kPipEDiff},
    {"EDiff", kEDiff},
    {"StartX", kStartX},
    {"StartY", kStartY},
    {"StartZ", kStartZ}
  };

  const std::unordered_map<KinematicTypes, std::string> KinematicParameterToStringMap = {
    {kTrueNeutrinoEnergy, "TrueNeutrinoEnergy"},
    {kyRec, "yRec"},
    {kRecoQ, "RecoQ"},
    {kRecoNeutrinoEnergy, "RecoNeutrinoEnergy"},
    {kIsFHC, "IsFHC"},
    {kyRec, "yRec"},
    {kOscChannel, "OscChannel"},
    {kMode, "Mode"},
    {kMuonMom, "MuonMom"},
    {kMuonEnergy, "MuonEnergy"},
    {kRecoMuonEnergy, "RecoMuonEnergy"},
    {kMuonTheta, "MuonTheta"},
    {kPipMom, "PipMom"},
    {kPipEnergy, "PipEnergy"},
    {kRecoPipEnergy, "RecoPipEnergy"},
    {kPipTheta, "PipTheta"},
    {kMuonEDiff, "MuonEDiff"},
    {kPipEDiff, "PipEDiff"},
    {kEDiff, "EDiff"},
    {kStartX, "StartX"},
    {kStartY, "StartY"},
    {kStartZ, "StartZ"}
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

  int _nuPDG;
  int _nuPDGunosc;

  double _PMom;
  double _PipMom;
  double _PimMom;
  double _Pi0Mom;
  double _NMom;
  double _MuMom;

  double _PTheta;
  double _PipTheta;
  double _PimTheta;
  double _Pi0Theta;
  double _NTheta;
  double _MuTheta;

  double _PStartX;
  double _PStartY;
  double _PStartZ;
  double _PEndX;
  double _PEndY;
  double _PEndZ;

  double _PipStartX;
  double _PipStartY;
  double _PipStartZ;
  double _PipEndX;
  double _PipEndY;
  double _PipEndZ;

  double _PimStartX;
  double _PimStartY;
  double _PimStartZ;
  double _PimEndX;
  double _PimEndY;
  double _PimEndZ;

  double _Pi0StartX;
  double _Pi0StartY;
  double _Pi0StartZ;
  double _Pi0EndX;
  double _Pi0EndY;
  double _Pi0EndZ;

  double _NStartX;
  double _NStartY;
  double _NStartZ;
  double _NEndX;
  double _NEndY;
  double _NEndZ;

  double _MuStartX;
  double _MuStartY;
  double _MuStartZ;
  double _MuEndX;
  double _MuEndY;
  double _MuEndZ;

  double _PMomReco;
  double _PipMomReco;
  double _PimMomReco;
  double _Pi0MomReco;
  double _NMomReco;
  double _MuMomReco;

  double _PThetaReco;
  double _PipThetaReco;
  double _PimThetaReco;
  double _Pi0ThetaReco;
  double _NThetaReco;
  double _MuThetaReco;

  double _eRecoP;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoN;
  double _eRecoMuon;

  double _LepNuAngle;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _eMuon;

  int _isCC;
  int nupdgUnosc;
  int nupdg;
  int _run;
  int _isND;
  int _isFHC;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;
  double _vtx_end_x;
  double _vtx_end_y;
  double _vtx_end_z;
  double _reco_vtx_x;
  double _reco_vtx_y;
  double _reco_vtx_z;
  double _reco_vtx_end_x;
  double _reco_vtx_end_y;
  double _reco_vtx_end_z;
  double _reco_px;
  double _reco_py;
  double _reco_pz;
  double _px;
  double _py;
  double _pz;
  double _LepTheta;
  double _Q2;
  int _reco_q;
  double _reco_pid;
  int _iscontained;

  double _E_diff;
  double _E_diff_Muon;
  double _E_diff_Pip;
  double _E_diff_Pim;
  double _E_diff_Pi0;
  double _E_diff_P;
  double _E_diff_N;
  // double* _mu_track_length_truth_true;
  // double* _mu_track_length_truth_reco;
  // double* _mu_track_length_reco;

  // configuration 
  bool IsELike;
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

  bool isNDCovSet = false;
  // The ND detector covariance matrix
  TMatrixD *NDCovMatrix;
  // The inverse ND detector covariance matrix
  double **NDInvertCovMatrix;


  std::vector<const double*> NDDetectorSystPointers;
  int nNDDetectorSystPointers;

  TH2* get2DParticleVarHist(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY);


};

#endif
