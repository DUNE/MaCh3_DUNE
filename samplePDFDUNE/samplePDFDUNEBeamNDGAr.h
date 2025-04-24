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
#include <random>
#include <limits>

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

    TH1* get1DParticleVarHist(std::string ProjectionVar_StrX, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* AxisX);
    TH2* get2DParticleVarHist(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY);
    
    enum KinematicTypes {kTrueNeutrinoEnergy, kRecoNeutrinoEnergy, kMode, kTrueXPos, kTrueYPos, kTrueZPos, kTrueRad, kNMuonsRecoOverTruth, kRecoLepEnergy, kTrueLepEnergy, kRecoXPos, kRecoYPos, kRecoZPos, kRecoRad, kLepPT, kLepPZ, kTrueQ0, kTrueQ3, kEvent_IsAccepted, kIsGoodCAFEvent, kParticle_Event, kParticle_Momentum, kParticle_TransverseMomentum, kParticle_BAngle, kParticle_IsAccepted, kParticle_PDG, kInFDV, kIsCC, kParticle_IsStoppedInTPC, kParticle_IsStoppedInECal, kParticle_IsStoppedInGap, kParticle_IsStoppedInEndGap, kParticle_NHits, kParticle_NTurns, kParticle_MomResMS, kParticle_MomResTrans};

  protected:
    void Init();
    int setupExperimentMC(int iSample);
    void setupFDMC(int iSample);

    void SetupWeightPointers();
    void SetupSplines();

    //DB functions which could be initialised to do something which is non-trivial
    //double CalcXsecWeightFunc(int iSample, int iEvent) {return 1.;}
    //void applyShifts(int iSample, int iEvent) {}

    const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
    const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
    const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);

    double ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent);
    double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
    double ReturnKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);

    std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);
    std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinematicParameter);

    void makePixelGrid(double pixel_spacing_cm);
    double FindNHits(double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature);
    double CalcBeta(double p_mag, double& bg, double& gamma);
    double GetMass(int partpdg);
    bool IsParticleAccepted(dunemc_base *duneobj, int i_sample, int i_event, int i_truepart, double pixel_spacing_cm, bool *isgoodcafparticle);

    bool IsParticleSelected(const int iSample, const int iEvent, const int iParticle);
    std::vector<struct dunemc_base> dunendgarmcSamples;

    TFile *_sampleFile;
    TTree *_data;
    TFile *_sampleFile_geant;
    TTree *_data_geant;
    TString _nutype;
    int _mode;

    double pot;
    int *nparticlesinsample;
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

    //Geant vectors
    std::vector<double> *_MCPStartX=0;
    std::vector<double> *_MCPStartY=0;
    std::vector<double> *_MCPStartZ=0;
    std::vector<double> *_MCPEndX=0;
    std::vector<double> *_MCPEndY=0;
    std::vector<double> *_MCPEndZ=0;
    std::vector<double> *_MCPStartPX=0;
    std::vector<double> *_MCPStartPY=0;
    std::vector<double> *_MCPStartPZ=0;
    std::vector<double> *_MCPEndPX=0;
    std::vector<double> *_MCPEndPY=0;
    std::vector<double> *_MCPEndPZ=0;
    std::vector<int> *_PDG = 0;
    std::vector<int> *_MCPTrkID=0;
    std::vector<int> *_SimHitTrkID=0;
    std::vector<double> *_SimHitEnergy=0;

    double pdgmass;
    //particle masses in GeV
    double m_chargedpi = 0.13957039;
    double m_pi0 = 0.1349768;
    double m_e = 0.00051099895;
    double m_mu = 0.1056583755;
    double m_p = 0.93827208816;
    double m_n = 0.9395654205;
    double m_chargedk = 0.493677;
    double m_k0 = 0.497648;
    double m_lambda = 1.115683;

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

    double X0 = 1193; //in cm From Federico's Kalman Filter Paper

    //pixel vars
    double pixelymin;
    double pixelymax;
    double pixelzmin;
    double pixelzmax;
    std::vector<double> yboundarypositions;
    std::vector<double> zboundarypositions;

    //configurable sample bools
    bool iscalo_reco; //NK Added so we can easily change what energy reconstruction we are using
    bool iselike;
    bool incl_geant; //NK - Added so we can use GArAnaTrees
    bool ecal_containment; //NK Do we count containment if the particle stops in the ECAL?

    //configurable sample cuts
    double muonscore_threshold; //NK Added so we can optimise muon threshold
    double protondEdxscore;
    double protontofscore;
    double recovertexradiusthreshold;
    double pionenergy_threshold; //NK Added so we can find pion energy threshold
    double B_field;
    double momentum_resolution_threshold;
    double pixel_spacing;
    double spatial_resolution;
    double adc_sampling_frequency;
    double drift_velocity;
    double pi0_reco_efficiency;  //efficiency for pi0 reco in ECAL 
    double gamma_reco_efficiency;  //efficiency for gamma reco in ECAL

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
      {"TrueQ0",kTrueQ0},
      {"TrueQ3",kTrueQ3},
      {"Event_IsAccepted",kEvent_IsAccepted},
      {"IsGoodCAFEvent",kIsGoodCAFEvent},
      {"Particle_Event",kParticle_Event},
      {"Particle_Momentum",kParticle_Momentum},
      {"Particle_TransverseMomentum",kParticle_TransverseMomentum},
      {"Particle_BAngle",kParticle_BAngle},
      {"Particle_IsAccepted",kParticle_IsAccepted},
      {"Particle_PDG",kParticle_PDG},
      {"InFDV",kInFDV},
      {"IsCC",kIsCC},
      {"Particle_IsStoppedInTPC",kParticle_IsStoppedInTPC},
      {"Particle_IsStoppedInECal",kParticle_IsStoppedInECal},
      {"Particle_IsStoppedInGap",kParticle_IsStoppedInGap},
      {"Particle_IsStoppedInEndGap",kParticle_IsStoppedInEndGap},
      {"Particle_NHits",kParticle_NHits},
      {"Particle_NTurns",kParticle_NTurns},
      {"Particle_MomResMS",kParticle_MomResMS},
      {"Particle_MomResTrans",kParticle_MomResTrans},
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
      {kTrueQ0,"TrueQ0"},
      {kTrueQ3,"TrueQ3"},
      {kEvent_IsAccepted,"Event_IsAccepted"},
      {kIsGoodCAFEvent,"IsGoodCAFEvent"},
      {kParticle_Event, "Particle_Event"},
      {kParticle_Momentum,"Particle_Momentum"},
      {kParticle_TransverseMomentum,"Particle_TransverseMomentum"},
      {kParticle_BAngle,"Particle_BAngle"},
      {kParticle_IsAccepted,"Particle_IsAccepted"},
      {kParticle_PDG,"Particle_PDG"},
      {kInFDV,"InFDV"},
      {kIsCC,"IsCC"},
      {kParticle_IsStoppedInTPC,"Particle_IsStoppedInTPC"},
      {kParticle_IsStoppedInECal,"Particle_IsStoppedInECal"},
      {kParticle_IsStoppedInGap,"Particle_IsStoppedInGap"},
      {kParticle_IsStoppedInEndGap,"Particle_IsStoppedInEndGap"},
      {kParticle_NHits,"Particle_NHits"},
      {kParticle_NTurns,"Particle_NTurns"},
      {kParticle_MomResMS,"Particle_MomResMS"},
      {kParticle_MomResTrans,"Particle_MomResTrans"},
    };
};

#endif
