#ifndef _samplePDFDUNEBeamFD_h_
#define _samplePDFDUNEBeamFD_h_

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

  enum KinematicTypes {kTrueNeutrinoEnergy,kRecoNeutrinoEnergy,kTrueXPos,kTrueYPos,kTrueZPos,kCVNNumu,kCVNNue,kM3Mode,kOscChannel,kIsFHC};

  //More robust getters to make plots in different variables, mode, osc channel, systematic weighting and with bin range 

  /// @brief Getter to make plots in different variables, mode, osc channel, systematic weighting and with bin range
  /// @param Var1 Variable kinematic type
  /// @param fModeToFill [default: -1] Mode to fill, -1 indicates all mdoes
  /// @param fSampleToFill [default: -1] Sample to fill, -1 indicates all samples
  /// @param WeightStyle [default: 0] If weight style==1 use *(MCSamples[i].osc_w_pointer[j]) * dunemcSamples[i].pot_s * dunemcSamples[i].norm_s * dunemcSamples[i].flux_w[j] instead
  /// @param Axis NOTE: UNUSED
  /// @return Plot of events binned in kinematic variable
  TH1D* get1DVarHist(KinematicTypes Var1, int fModeToFill=-1, int fSampleToFill=-1, int WeightStyle=0, TAxis* Axis=0);

  /// @brief Getter to make plots in different variables, mode, osc channel, systematic weighting and with bin range
  /// @param Var1 Variable kinematic type
  /// @param Selection Selection to apply to hist
  /// @param WeightStyle [default: 0] If weight style==1 use *(MCSamples[i].osc_w_pointer[j]) * dunemcSamples[i].pot_s * dunemcSamples[i].norm_s * dunemcSamples[i].flux_w[j] instead
  /// @param Axis NOTE: UNUSED
  /// @return Plot of events binned in kinematic variable
  TH1D* get1DVarHist(KinematicTypes Var1, std::vector< std::vector<double> > Selection, int WeightStyle=0, TAxis* Axis=0);
    
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
  
  /// @todo extract this completely to core
  ///@brief Setup our spline file, this calls InitialseSplineObject() under the hood
  void SetupSplines();
  
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event 
  double ReturnKinematicParameter (double KinematicVariable, int iSample, int iEvent);

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
  inline int ReturnKinematicParameterFromString(std::string KinematicParameterStr);
  inline std::string ReturnStringFromKinematicParameter(int KinematicParameterStr);
  
  //DB functions which could be initialised to do something which is non-trivial
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iSample sample ID
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}

  /// @brief Apply kinematic shifts
  /// @param iSample Sample Number
  /// @param iEvent Event number
  void applyShifts(int iSample, int iEvent);

  // dunemc
  /// DUNE MC sampels
  std::vector<struct dunemc_base> dunemcSamples;

  /// Value of POT used for sample
  double pot;

  /// File containing sample objects
  TFile *_sampleFile;

  /// TTree containing sample Data
  TTree *_data;

  TString _nutype;

  int _mode;

  //Reco Variables
  double _erec;
  double _erec_nue;
  double _erec_had;
  double _erec_had_nue;
  double _erec_lep;
  double _erec_lep_nue;

  double _eRecoP;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoN;

  double _cvnnumu;
  double _cvnnue;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;

  //Truth Variables
  double _ev;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _NuMomX;
  double _NuMomY;
  double _NuMomZ;
  double _LepMomX;
  double _LepMomY;
  double _LepMomZ;
  double _LepNuAngle;
  double _BeRPA_cvwgt;
  double _maccres_cvwgt;
  double _nunpcc1_cvwgt;
  int _isCC;
  int _nuPDGunosc;
  int _nuPDG;
  int _run;
  int _isFHC;
  double _LepTheta;
  double _Q2;

  bool iselike;

  //Positions of FD Detector systematics
  double tot_escale_fd_pos;
  double tot_escale_sqrt_fd_pos;
  double tot_escale_invsqrt_fd_pos;
  double had_escale_fd_pos;
  double had_escale_sqrt_fd_pos;
  double had_escale_invsqrt_fd_pos;
  double mu_escale_fd_pos;
  double mu_escale_sqrt_fd_pos;
  double mu_escale_invsqrt_fd_pos;
  double n_escale_fd_pos;
  double n_escale_sqrt_fd_pos;
  double n_escale_invsqrt_fd_pos;
  double em_escale_fd_pos;
  double em_escale_sqrt_fd_pos;
  double em_escale_invsqrt_fd_pos;
  double had_res_fd_pos;
  double mu_res_fd_pos;
  double n_res_fd_pos;
  double em_res_fd_pos;
  double cvn_numu_fd_pos;
  double cvn_nue_fd_pos;

  /// FD Detector Systematics
  std::vector<const double*> FDDetectorSystPointers;

  /// Number of FD Detector Systematics
  int nFDDetectorSystPointers;
};

#endif
