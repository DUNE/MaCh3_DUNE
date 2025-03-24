#ifndef _samplePDFDUNEBeamNDGAr_h_
#define _samplePDFDUNEBeamNDGAr_h_

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

class samplePDFDUNEBeamNDGAr : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBeamNDGAr(std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBeamNDGAr();

  enum KinematicTypes {kTrueNeutrinoEnergy, kRecoNeutrinoEnergy, kMode, kTrueXPos, kTrueYPos, kTrueZPos, kTrueRad, kNMuonsRecoOverTruth, kRecoLepEnergy, kTrueLepEnergy, kRecoXPos, kRecoYPos, kRecoZPos, kRecoRad, kLepPT, kLepPZ, kPionMultiplicity, kNRecoParticles, kInFDV, kTrueMinusRecoEnergyRatio, kTrueMinusRecoEnergy, kNTrueMuons, kNRecoMuons, kOscChannel, kTrueQ0, kTrueQ3};
  
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
  inline std::string ReturnStringFromKinematicParameter(int KinematicParameter);
  
	void makePixelGrid(float pixel_spacing_cm);
	double FindNHits(float pixel_spacing_cm, float centre_circle_y, float centre_circle_z, double rad_curvature);
  double CalcBeta(double p_mag, double& bg, double& gamma);
	double CalcDeDx(double beta, double bg, double gamma);

	// dunendmc
  std::vector<struct dunemc_base> dunendgarmcSamples;

  TFile *_sampleFile;
  TTree *_data;
	TFile *_sampleFile_geant;
	TTree *_data_geant;
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

	double pdgmass;
  //particle masses in GeV
  double m_chargedpi = 0.13957039;
  double m_pi0 = 0.1349768;
  double m_e = 0.00051099895;
  double m_mu = 0.1056583755;
  double m_p = 0.93827208816;
  double m_n = 0.9395654205;
  double m_chargedk = 0.493677;

	//TPC dimensions
	double TPCFidLength;
  double TPCFidRadius;
  double TPCInstrumentedLength;
  double TPCInstrumentedRadius;
  double ECALInnerRadius;
  double ECALOuterRadius;
  double ECALEndCapStart;
  double ECALEndCapEnd;

	double TPC_centre_x =0.;
  double TPC_centre_y = -150.;
  double TPC_centre_z = 1486.;

	double K_const = 0.307075; //4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol)
  double sternheimer_A = 0.1956;
  double sternheimer_K = 3.0000;
  double sternheimer_X0 = 0.2000;
  double sternheimer_X1 = 3.0000;
  double sternheimer_Cbar = 5.2146;
  double excitationenergy = 188.0; //excitation energy for electrons in argon gas in eV
  double density = 0.0167; //in g/cm^3
  double X0 = 1193; //in cm From Federico's Kalman Filter Paper

	//pixel vars
	double pixelymin;
  double pixelymax;
  double pixelzmin;
  double pixelzmax;
  std::vector<double> yboundarypositions;
  std::vector<double> zboundarypositions;

  bool iscalo_reco; //NK Added so we can easily change what energy reconstruction we are using
  bool iselike;
	bool incl_geant; //NK - Added so we can use GArAnaTrees
	bool ecal_containment; //NK Do we count containment if the particle stops in the ECAL?

	float muonscore_threshold; //NK Added so we can optimise muon threshold
	float protondEdxscore;
  float protontofscore;
  float recovertexradiusthreshold;
  float pionenergy_threshold; //NK Added so we can find pion energy threshold
  float B_field;
  float momentum_resolution_threshold;
  float pixel_spacing;
  float spatial_resolution;
  float adc_sampling_frequency;
  float drift_velocity;
//  float hits_per_mm;
  float pi0_reco_efficiency;  //efficiency for pi0 reco in ECAL 
  float gamma_reco_efficiency;  //efficiency for gamma reco in ECAL

  caf::StandardRecord* sr = new caf::StandardRecord();

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"RecoNeutrinoEnergy",kRecoNeutrinoEnergy},
	{"Mode",kMode},
    {"TrueXPos",kTrueXPos},
    {"TrueYPos",kTrueYPos},
    {"TrueZPos",kTrueZPos},
    {"TrueRad",kTrueRad},
    {"NMuonsRecoOverTruth",kNMuonsRecoOverTruth},
    {"RecoLepEnergy",kRecoLepEnergy},
    {"TrueLepEnergy",kTrueLepEnergy},
	{"RecoXPos",kRecoXPos},
	{"RecoYPos",kRecoYPos},
	{"RecoZPos",kRecoZPos},
	{"RecoRad",kRecoRad},
	{"LepPT",kLepPT},
	{"LepPZ",kLepPZ},
	{"OscillationChannel",kOscChannel},
	{"TrueQ0",kTrueQ0},
	{"TrueQ3",kTrueQ3}
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kRecoNeutrinoEnergy,"RecoNeutrinoEnergy"},
	{kMode,"Mode"},
    {kTrueXPos,"TrueXPos"},
    {kTrueYPos,"TrueYPos"},
	{kTrueZPos,"TrueZPos"},
	{kTrueRad,"TrueRad"},
    {kNMuonsRecoOverTruth,"NMuonsRecoOverTruth"},
    {kRecoLepEnergy,"RecoLepEnergy"},
    {kTrueLepEnergy,"TrueLepEnergy"},
    {kRecoXPos,"RecoXPos"},
    {kRecoYPos,"RecoYPos"},
    {kRecoZPos,"RecoZPos"},
    {kRecoRad,"RecoRad"},
    {kLepPT,"LepPT"},
    {kLepPZ,"LepPZ"},
	{kOscChannel,"OscillationChannel"},
	{kTrueQ0,"TrueQ0"},
	{kTrueQ3,"TrueQ3"},
  };
};

#endif
