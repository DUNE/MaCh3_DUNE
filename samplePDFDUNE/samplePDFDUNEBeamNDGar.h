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
#include "StandardRecord.h"

#include "StructsDUNE.h"

/// @brief Base class for handling beam ND GAR samples
class samplePDFDUNEBeamNDGar : virtual public samplePDFFDBase
{
public:
  /// @brief Constructor
  /// @param mc_version Configuration file
  /// @param xsec_cov cross-section covariance matrix
  samplePDFDUNEBeamNDGar(std::string mc_version, covarianceXsec* xsec_cov);

  /// @brief destructor
  ~samplePDFDUNEBeamNDGar();

  /// @brief Enum to identify kinematics
  enum KinematicTypes {kTrueNeutrinoEnergy, kRecoNeutrinoEnergy, kTrueXPos, kTrueYPos, kTrueZPos, kTrueRad, kNMuonsRecoOverTruth, kRecoLepEnergy, kTrueLepEnergy, kRecoXPos, kRecoYPos, kRecoZPos, kRecoRad, kLepPT, kLepPZ, kPionMultiplicity, kNRecoParticles, kInFDV, kTrueMinusRecoEnergyRatio, kTrueMinusRecoEnergy, kNTrueMuons, kNRecoMuons};
  
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

  /// @brief Sets up splines 
  void SetupSplines();
  
  /// @brief NOT IMPLEMENTED: Dunder method to calculate xsec weights
  /// @param iSample sample ID
  /// @param iEvent Event number
  double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}

  /// @brief NOT IMPLEMENTED: Apply kinematic shifts
  /// @param iSample Sample Number
  /// @param iEvent Event number
  void applyShifts(int iSample, int iEvent) {}

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinPar Kinematic parameter enum val
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Pointer to KinPar for a given event
  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinetmatic parameter corresponding for a given event
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicVariable Kinematic parameter ID as double (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinetmatic parameter corresponding for a given event
  double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
  
  /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
  /// @param KinematicParameter Kinematic parameter name as string (gets cast -> int)
  /// @param iSample Sample ID
  /// @param iEvent Event ID
  /// @return Value of kinematic parameter corresponding for a given event
  double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
  
  /// @brief Gets binning for a given parameter
  /// @param KinematicParameterStr Parameter name
  /// @return Vector containing parameter bins
  std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);

  /// @brief Gets binning for a given parameter
  /// @param KinPar Parameter ID
  /// @return Vector containing parameter bins
  std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinematicParameter);

  /// @brief Get kinematic parameter ID from string name
  /// @param KinematicStr 
  /// @return Parameter ID
  int ReturnKinematicParameterFromString(std::string KinematicParameterStr);

  /// @brief Gets name of kinematic parmaeter
  /// @param KinPar Parameter ID
  /// @return Name of parameter
  std::string ReturnStringFromKinematicParameter(int KinematicParameter);
  
  // dunendmc
  /// Array filled with MC samples for each oscillation channel
  std::vector<struct dunemc_base> dunendgarmcSamples;

  /// File containing sample objects
  TFile *_sampleFile;
  /// TTree containing sample Data
  TTree *_data;
  /// Value of POT used for sample
  double pot;


  TString _nutype;
  int _mode;


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
