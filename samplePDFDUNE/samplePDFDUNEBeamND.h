#ifndef _samplePDFDUNEBeamND_h_
#define _samplePDFDUNEBeamND_h_

#include <cstddef>
#include <iostream>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraph2DErrors.h>
#include <manager/MaCh3Modes.h>
#include <vector>
#include <omp.h>
#include <list>

#include "splines/splinesDUNE.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"
#include "samplePDF/Structs.h"

class samplePDFDUNEBeamND : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamND(std::string mc_version, covarianceXsec* xsec_cov, covarianceOsc* osc_cov);
  ~samplePDFDUNEBeamND();

  enum KinematicTypes {kTrueNeutrinoEnergy, kRecoQ, kRecoNeutrinoEnergy, kIsFHC, kRecoY};
  
 protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();

  // === HH: Functional parameters ===
  enum FuncParEnum {kTotalEScaleND, kDebugNothing, kDebugShift};
  void RegisterFunctionalParameters();
  void SetupFunctionalParameters();

  // Map the name of the functional parameter to funcpar enum
  std::unordered_map<std::string, FuncParEnum> funcParsNamesMap;

  // Map the funcpar enum to pointer of FuncPars struct
  std::unordered_map<FuncParEnum, FuncPars*> funcParsMap;

  std::unordered_map<FuncParEnum, FuncParFuncType> funcParsFuncMap;

  // Map the funcpar enum to the value of the functional parameter
  std::unordered_map<FuncParEnum, double*> funcParsValMap;

  // A grid of vectors of enums for each sample and event
  std::vector<std::vector<std::vector<FuncParEnum>>> funcParsGrid;

  void TotalEScale(const double * par, std::size_t iSample, std::size_t iEvent);
  void DebugShift(const double * par, std::size_t iSample, std::size_t iEvent);
  // =================================
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  int ReturnKinematicParameterFromString(std::string KinematicParameterStr);
  std::string ReturnStringFromKinematicParameter(int KinematicParameter);
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
  void applyShifts(int iSample, int iEvent);

  std::vector<struct dunemc_base> dunendmcSamples;

  TFile *_sampleFile;
  TTree *_data;
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
  bool IsRHC;

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

  int nNDDetectorSystPointers = 0;
};

#endif
