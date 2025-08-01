#ifndef _SampleHandlerBeamNDGAr_h_
#define _SampleHandlerBeamNDGAr_h_

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerFD.h"
#include "Samples/StructsDUNE.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

class SampleHandlerBeamNDGAr : virtual public SampleHandlerFD
{
  public:
    SampleHandlerBeamNDGAr(std::string mc_version, ParameterHandlerGeneric* xsec_cov);
    ~SampleHandlerBeamNDGAr();

    TH1* Get1DParticleVarHist(std::string ProjectionVar_StrX, std::vector< KinematicCut > SelectionVec, int WeightStyle, TAxis* AxisX);
    TH2* Get2DParticleVarHist(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, std::vector< KinematicCut > SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY);
    
    enum KinematicTypes {kTrueNeutrinoEnergy, kRecoNeutrinoEnergy, kMode, kTrueXPos, kTrueYPos, kTrueZPos, kTrueRad, kNMuonsRecoOverTruth, kRecoLepEnergy, kTrueLepEnergy, kRecoXPos, kRecoYPos, kRecoZPos, kRecoRad, kLepPT, kLepPZ, kLepP, kLepBAngle, kLepTheta, kLepPhi, kTrueQ0, kTrueQ3, kEvent_IsAccepted, kIsGoodCAFEvent, kParticle_Event, kParticle_Momentum, kParticle_EndMomentum, kParticle_TransverseMomentum, kParticle_BAngle, kParticle_BeamAngle, kParticle_IsAccepted, kParticle_IsContained, kParticle_IsDecayed, kParticle_PDG, kInFDV, kIsCC, kParticle_IsStoppedInTPC, kParticle_IsStoppedInECal, kParticle_IsStoppedInGap, kParticle_IsStoppedInEndGap, kParticle_IsStoppedInBarrelGap, kParticle_IsEscaped, kParticle_NTurns, kParticle_NHits, kParticle_TrackLengthYZ, kParticle_MomResMS, kParticle_MomResYZ, kParticle_MomResX, kParticle_StartR2, kParticle_EndR, kParticle_EndX, kParticle_StartX, kParticle_EDepCrit, kParticle_EscSecEnergyFrac};

  protected:
    //Functions required by core
    void Init();
    int SetupExperimentMC(int iSample);
    void SetupFDMC(int iSample);

    void SetupWeightPointers();
    void SetupSplines();

    void RegisterFunctionalParameters() override {};

    const double* GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent);
    const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
    double ReturnKinematicParameter(int KinematicVariable, int iSample, int iEvent);
    double ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent);
    std::vector<double> ReturnKinematicParameterBinning(std::string KinematicParameter);

    //NDGAr-specific functions
    const double* GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
    double ReturnKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent);
    std::vector<double> ReturnKinematicParameterBinning(KinematicTypes KinematicParameter);
    
    double FindNHits(double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature, double theta_start, double theta_spanned, int charge);
    bool isCoordOnTrack(int charge, double ycoord, double zcoord, double centre_circle_y, double centre_circle_z, double theta_start, double theta_spanned);
    double CalcBeta(double p_mag, double& bg, double& gamma, double pdgmass);
    int GetChargeFromPDG(int pdg);
    bool IsParticleAccepted(dunemc_base *duneobj, int i_sample, int i_event, size_t i_anapart, double pixel_spacing_cm, 
         std::unordered_map<int,std::vector<double>>& primID_to_EDepCrit, std::unordered_map<int, std::vector<int>>& prim_to_sec_ID, std::unordered_map<int, size_t>& ID_to_index);

    bool IsParticleSelected(const int iSample, const int iEvent, const int iParticle);
    std::vector<struct dunemc_base> dunendgarmcSamples;

    TFile *_sampleFile;
    TTree *_data;
    TFile *_sampleFile_geant;
    TTree *_data_geant;
    TString _nutype;
    //int _mode;

    double pot;
    int *nparticlesinsample;
    double _BeRPA_cvwgt = 1;

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
    std::vector<std::string> *_MCPProc=0;
    std::vector<std::string> *_MCPEndProc=0;
    std::vector<int> *_MotherTrkID=0;
    std::vector<int> *_SimHitTrkID=0;
    std::vector<int> *_SimHitLayer=0;
    std::vector<double> *_SimHitEnergy=0;
    std::vector<double> *_SimHitX=0;
    std::vector<double> *_SimHitY=0;
    std::vector<double> *_SimHitZ=0;

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
    double BeamDirection[3] = {0.,-0.101,0.995};

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
    std::ofstream deposit_outputs;

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
      {"LepTheta",kLepTheta},
      {"LepPhi",kLepPhi},
      {"LepP",kLepP},
      {"LepBAngle",kLepBAngle},
      {"TrueQ0",kTrueQ0},
      {"TrueQ3",kTrueQ3},
      {"Event_IsAccepted",kEvent_IsAccepted},
      {"IsGoodCAFEvent",kIsGoodCAFEvent},
      {"Particle_Event",kParticle_Event},
      {"Particle_Momentum",kParticle_Momentum},
      {"Particle_EndMomentum",kParticle_EndMomentum},
      {"Particle_TransverseMomentum",kParticle_TransverseMomentum},
      {"Particle_BAngle",kParticle_BAngle},
      {"Particle_BeamAngle",kParticle_BeamAngle},
      {"Particle_IsAccepted",kParticle_IsAccepted},
      {"Particle_IsContained",kParticle_IsContained},
      {"Particle_IsDecayed",kParticle_IsDecayed},
      {"Particle_PDG",kParticle_PDG},
      {"InFDV",kInFDV},
      {"IsCC",kIsCC},
      {"Particle_IsStoppedInTPC",kParticle_IsStoppedInTPC},
      {"Particle_IsStoppedInECal",kParticle_IsStoppedInECal},
      {"Particle_IsStoppedInGap",kParticle_IsStoppedInGap},
      {"Particle_IsStoppedInEndGap",kParticle_IsStoppedInEndGap},
      {"Particle_IsStoppedInBarrelGap",kParticle_IsStoppedInBarrelGap},
      {"Particle_IsEscaped",kParticle_IsEscaped},
      {"Particle_NTurns",kParticle_NTurns},
      {"Particle_NHits",kParticle_NHits},
      {"Particle_TrackLengthYZ",kParticle_TrackLengthYZ},
      {"Particle_MomResMS",kParticle_MomResMS},
      {"Particle_MomResYZ",kParticle_MomResYZ},
      {"Particle_MomResX",kParticle_MomResX},
      {"Particle_StartR2",kParticle_StartR2},
      {"Particle_EndR",kParticle_EndR},
      {"Particle_EndX",kParticle_EndX},
      {"Particle_StartX",kParticle_StartX},
      {"Particle_EDepCrit",kParticle_EDepCrit},
      {"Particle_EscSecEnergyFrac",kParticle_EscSecEnergyFrac},
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
      {kLepTheta,"LepTheta"},
      {kLepPhi,"LepPhi"},
      {kLepBAngle,"LepBAngle"},
      {kLepP,"LepP"},
      {kTrueQ0,"TrueQ0"},
      {kTrueQ3,"TrueQ3"},
      {kEvent_IsAccepted,"Event_IsAccepted"},
      {kIsGoodCAFEvent,"IsGoodCAFEvent"},
      {kParticle_Event, "Particle_Event"},
      {kParticle_Momentum,"Particle_Momentum"},
      {kParticle_EndMomentum,"Particle_EndMomentum"},
      {kParticle_TransverseMomentum,"Particle_TransverseMomentum"},
      {kParticle_BAngle,"Particle_BAngle"},
      {kParticle_BeamAngle,"Particle_BeamAngle"},
      {kParticle_IsAccepted,"Particle_IsAccepted"},
      {kParticle_IsContained,"Particle_IsContained"},
      {kParticle_IsDecayed,"Particle_IsDecayed"},
      {kParticle_PDG,"Particle_PDG"},
      {kInFDV,"InFDV"},
      {kIsCC,"IsCC"},
      {kParticle_IsStoppedInTPC,"Particle_IsStoppedInTPC"},
      {kParticle_IsStoppedInECal,"Particle_IsStoppedInECal"},
      {kParticle_IsStoppedInGap,"Particle_IsStoppedInGap"},
      {kParticle_IsStoppedInEndGap,"Particle_IsStoppedInEndGap"},
      {kParticle_IsStoppedInBarrelGap,"Particle_IsStoppedInBarrelGap"},
      {kParticle_IsEscaped,"Particle_IsEscaped"},
      {kParticle_NTurns,"Particle_NTurns"},
      {kParticle_NHits,"Particle_NHits"},
      {kParticle_TrackLengthYZ,"Particle_TrackLengthYZ"},
      {kParticle_MomResMS,"Particle_MomResMS"},
      {kParticle_MomResYZ,"Particle_MomResYZ"},
      {kParticle_MomResX,"Particle_MomResX"},
      {kParticle_StartR2,"Particle_StartR2"},
      {kParticle_EndR,"Particle_EndR"},
      {kParticle_EndX,"Particle_EndX"},
      {kParticle_StartX,"Particle_StartX"},
      {kParticle_EDepCrit,"Particle_EDepCrit"},
      {kParticle_EscSecEnergyFrac,"Particle_EscSecEnergyFrac"},
    };
};

#endif
