#pragma once

#include <vector>

#include "covariance/covarianceOsc.h"
#include "covariance/covarianceXsec.h"

#include "samplePDF/samplePDFFDBase.h"

#include "StructsDUNE.h"

class samplePDFDUNEBeamFD : virtual public samplePDFFDBase {
public:
  samplePDFDUNEBeamFD(std::string mc_version, covarianceXsec *xsec_cov);
  ~samplePDFDUNEBeamFD();

  double getPurity_numu();
  double getEfficiency_numu();
  double getCut();

  double getPurity_nue();
  double getEfficiency_nue();

  enum KinematicTypes {
    kTrueNeutrinoEnergy,
    kRecoNeutrinoEnergy,
    kTrueXPos,
    kTrueYPos,
    kTrueZPos,
    kCVNNumu,
    kCVNNue,
    kM3Mode,
    kOscChannel,
    kNumu_efficiency,
    kNue_efficiency,
    kselected_numuCCevent_energy,
    kselected_numuCCevent_vertexpos_x,
    kselected_numuCCevent_vertexpos_y,
    kselected_numuCCevent_vertexpos_z,
    kselected_nueCCevent_energy,
    kselected_nueCCevent_vertexpos_x,
    kselected_nueCCevent_vertexpos_y,
    kselected_nueCCevent_vertexpos_z,
    kismc_numu,
    kismc_nue,
    kismc_nutau,
    kisCC,
    knuflavour,
    krecopdg

  };

  // More robust getters to make plots in different variables, mode, osc
  // channel, systematic weighting and with bin range
  TH1D *get1DVarHist(KinematicTypes Var1, int fModeToFill = -1,
                     int fSampleToFill = -1, int WeightStyle = 0,
                     TAxis *Axis = 0);
  TH1D *get1DVarHist(KinematicTypes Var1,
                     std::vector<std::vector<double>> Selection,
                     int WeightStyle = 0, TAxis *Axis = 0);

protected:
  void Init();
  int setupExperimentMC(int iSample);
  void setupFDMC(int iSample);

  void SetupWeightPointers();

  /// @todo extract this completely to core
  ///@brief Setup our spline file, this calls InitialseSplineObject() under the
  /// hood
  void SetupSplines();

  double CalculatePOT();

  double ReturnKinematicParameter(double KinematicVariable, int iSample,
                                  int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample,
                                  int iEvent);

  const double *GetPointerToKinematicParameter(std::string KinematicParameter,
                                               int iSample, int iEvent);
  const double *GetPointerToKinematicParameter(double KinematicVariable,
                                               int iSample, int iEvent);

  std::vector<double>
  ReturnKinematicParameterBinning(std::string KinematicParameter);
  inline int
  ReturnKinematicParameterFromString(std::string KinematicParameterStr);
  inline std::string
  ReturnStringFromKinematicParameter(int KinematicParameterStr);

  // DB functions which could be initialised to do something which is
  // non-trivial
  double CalcXsecWeightFunc(int iSample, int iEvent) { return 1.; }
  void applyShifts(int iSample, int iEvent);

  // dunemc
  std::vector<dunemc_base> dunemcSamples;

  double gen_pot;
  double numu_cut;
  double nue_cut;
};
