#ifndef _samplePDFDUNEAtm_h_
#define _samplePDFDUNEAtm_h_

#include "splines/splinesDUNE.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"
#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"

class samplePDFDUNEAtm : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEAtm(std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEAtm();

  enum KinematicTypes{kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueCosZ,kRecoCosZ};
  
 protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();
  void SetupSplines();
  
  //DB functions which could be initialised to do something which is non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
  void applyShifts(int iSample, int iEvent) {}
  
  double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
  double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
  
  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
  
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
  inline int ReturnKinematicParameterFromString(std::string KinematicStr);
  inline std::string ReturnStringFromKinematicParameter(int KinematicVariable);

  std::vector<struct dunemc_base> dunemcSamples;
  bool IsELike;
};

#endif
