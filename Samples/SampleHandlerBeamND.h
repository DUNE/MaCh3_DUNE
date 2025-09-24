#ifndef _SampleHandlerBeamND_h_
#define _SampleHandlerBeamND_h_

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerFD.h"

#include "StructsDUNE.h"

class SampleHandlerBeamND : virtual public SampleHandlerFD
{
public:
  SampleHandlerBeamND(std::string mc_version, ParameterHandlerGeneric* xsec_cov, TMatrixD* nd_cov) ;
  ~SampleHandlerBeamND();

  enum KinematicTypes {kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kyRec,kOscChannel,kMode,kIsFHC};
  
 protected:
  void Init();
  int SetupExperimentMC();
  void SetupFDMC();

  void SetupWeightPointers();
  void SetupSplines();

  void RegisterFunctionalParameters() override {};
  
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent);
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent);

  double ReturnKinematicParameter(int KinematicVariable, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent);

  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iEvent) {return 1.; (void)iEvent;}

  void setNDCovMatrix();
  double GetLikelihood() override;

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
  double pot_s;
  double norm_s;

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

  bool isNDCovSet = false;
  // The ND detector covariance matrix
  TMatrixD *NDCovMatrix;
  // The inverse ND detector covariance matrix
  double **NDInvertCovMatrix;


  std::vector<const double*> NDDetectorSystPointers;
  int nNDDetectorSystPointers;
  std::unordered_map<std::string, std::vector<double>> norm_map;

  /// @brief Cleanup memory
  void CleanMemoryBeforeFit() override {};
};

#endif
