#ifndef _samplePDFDUNEBeamNDGar_h_
#define _samplePDFDUNEBeamNDGar_h_

#include <iostream>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraph2DErrors.h>
#include <vector>
#include <omp.h>
#include <list>

#include "splines/splinesDUNE.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"
#include "samplePDF/samplePDFFDBase.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "StructsDUNE.h"

class samplePDFDUNEBeamNDGar : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamNDGar(std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBeamNDGar();

  enum KinematicTypes {kTrueNeutrinoEnergy, kRecoNeutrinoEnergy, kTrueXPos, kTrueYPos, kTrueZPos, kTrueRad, kNMuonsRecoOverTruth, kRecoLepEnergy, kTrueLepEnergy, kRecoXPos, kRecoYPos, kRecoZPos, kRecoRad, kLepPT, kLepPZ, kPionMultiplicity, kNRecoParticles, kInFDV, kTrueMinusRecoEnergyRatio, kTrueMinusRecoEnergy, kNTrueMuons, kNRecoMuons};
  
 protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
  void applyShifts(int iSample, int iEvent) {}

  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinematicParameter);
  int ReturnKinematicParameterFromString(std::string KinematicParameterStr);
  std::string ReturnStringFromKinematicParameter(int KinematicParameter);
  
  // dunendmc
  std::vector<struct dunemc_base> dunendgarmcSamples;

  TFile *_sampleFile;
  TTree *_data;
  TString _nutype;
  int _mode;

  double pot;

  // dunendgarmc Variables
  double _ev;
  double _erec;
  double _erec_nue;
  double _elep_reco;
  double _LepNuAngle;
  int _reco_numu;
  int _reco_nue;
  double _BeRPA_cvwgt = 1;
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

  bool iscalo_reco; //NK Added so we can easily change what energy reconstruction we are using
  float muonscore_threshold; //NK Added so we can optimise muon threshold

  caf::StandardRecord* sr = new caf::StandardRecord();
};

#endif
