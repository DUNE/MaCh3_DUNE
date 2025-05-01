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

  TH2* get2DParticleVarHist(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY);

  enum KinematicTypes {
    kTrueNeutrinoEnergy, 
    kRecoQ, 
    kRecoNeutrinoEnergy, 
    kIsFHC,
    kyRec, 
    kOscChannel, 
    kMode, 
    kParticle_MuonMom, 
    kParticle_MuonEnergy, 
    kRecoMuonEnergy, 
    kParticle_MuonTheta, 
    kParticle_PipMom,
    kParticle_PipEnergy,
    kRecoPipEnergy, 
    kParticle_PipTheta, 
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
    {"MuonMom", kParticle_MuonMom},
    {"MuonEnergy", kParticle_MuonEnergy},
    {"MuonTheta", kParticle_MuonTheta},
    {"RecoMuonEnergy", kRecoMuonEnergy},
    {"PipMom", kParticle_PipMom},
    {"PipEnergy", kParticle_PipEnergy},
    {"RecoPipEnergy", kRecoPipEnergy},
    {"PipTheta", kParticle_PipTheta},
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
    {kParticle_MuonMom, "Particle_MuonMom"},
    {kParticle_MuonEnergy, "Particle_MuonEnergy"},
    {kRecoMuonEnergy, "RecoMuonEnergy"},
    {kParticle_MuonTheta, "Particle_MuonTheta"},
    {kParticle_PipMom, "Particle_PipMom"},
    {kParticle_PipEnergy, "Particle_PipEnergy"},
    {kRecoPipEnergy, "RecoPipEnergy"},
    {kParticle_PipTheta, "Particle_PipTheta"},
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
  int *nparticlesinsample;
  
  bool isNDCovSet = false;
  // The ND detector covariance matrix
  TMatrixD *NDCovMatrix;
  // The inverse ND detector covariance matrix
  double **NDInvertCovMatrix;


  std::vector<const double*> NDDetectorSystPointers;
  int nNDDetectorSystPointers;

};

#endif
