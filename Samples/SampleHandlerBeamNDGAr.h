#ifndef _SampleHandlerBeamNDGAr_h_
#define _SampleHandlerBeamNDGAr_h_

#include "Samples/SampleHandlerFD.h"
#include "Samples/StructsDUNE.h"
#include "Splines/BinnedSplineHandlerDUNE.h"

class SampleHandlerBeamNDGAr : virtual public SampleHandlerFD
{
public:
  SampleHandlerBeamNDGAr(std::string mc_version, ParameterHandlerGeneric* xsec_cov);
  ~SampleHandlerBeamNDGAr();

  enum KinematicTypes {kTrueNeutrinoEnergy, kMode, kOscChannel, kTrueXPos, kTrueYPos, kTrueZPos, kTrueRad, kTrueLepEnergy,
    kLepPT, kLepPZ, kLepP, kLepBAngle, kLepTheta, kLepPhi, kTrueQ0, kTrueQ3, kEvent_IsAccepted,  kInFDV, kIsCC, kEPi0, kNPi0,
    kLepTrackLengthYZ};
  
  enum KinematicVecs {kParticle_EVis, kParticle_Momentum, kParticle_EndMomentum, kParticle_TransverseMomentum, 
    kParticle_BAngle, kParticle_BeamAngle, kParticle_IsAccepted, kParticle_IsCurvatureResolved, kParticle_IsDecayed, kParticle_PDG,
    kParticle_IsStoppedInTPC, kParticle_IsStoppedInECal, kParticle_IsStoppedInBarrel, kParticle_IsStoppedInEndCap, kParticle_IsStoppedInGap, 
    kParticle_IsStoppedInEndGap, kParticle_IsStoppedInBarrelGap, kParticle_IsEscaped, kParticle_NTurns, kParticle_NHits,
    kParticle_TrackLengthYZ, kParticle_MomResMS, kParticle_MomResYZ, kParticle_MomResX, kParticle_StartR2, kParticle_EndR, 
    kParticle_EndDepth, kParticle_EndX, kParticle_EndY, kParticle_EndZ, kParticle_StartX, kParticle_EDepCrit};

protected:
  //Functions required by core
  void Init() override;
  int SetupExperimentMC() override;
  void SetupFDMC() override;
  void SetupWeightPointers() override;
  void SetupSplines() override;

  void CleanMemoryBeforeFit() override;
  void RegisterFunctionalParameters() override {};

  const double* GetPointerToKinematicParameter(KinematicTypes KinPar, size_t iEvent);
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) override;
  const double* GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) override;

  double ReturnKinematicParameter(KinematicTypes KinPar, size_t iEvent); // Extra function to deal with non-doubles
  double ReturnKinematicParameter(int KinematicVariable, int iEvent) override;
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent) override;

  std::vector<double> ReturnKinematicVector(KinematicVecs KinVec, size_t iEvent);
  std::vector<double> ReturnKinematicVector(int KinematicVector, int iEvent) override;
  std::vector<double> ReturnKinematicVector(std::string KinematicVector, int iEvent) override;

  std::vector<dunemc_beamndgar> dunendgarmcFitting;
  std::vector<dunemc_plotting> dunendgarmcPlotting;

  //NDGAr-specific functions
  double FindNHits(double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature, double theta_start, double theta_spanned, int charge);
  bool isCoordOnTrack(int charge, double ycoord, double zcoord, double centre_circle_y, double centre_circle_z, double theta_start, double theta_spanned);
  double CalcBeta(double p_mag, double& bg, double& gamma, double pdgmass);
  int GetChargeFromPDG(int pdg);
  bool IsResolvedFromCurvature(dunemc_plotting& plotting_vars, size_t i_anapart, double pixel_spacing_cm);
  double GetCalDepth(double x, double y, double z);
  double DepthToLayer(double depth, double r);
  double CalcEDepCal(int motherID, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, std::vector<double>>& ID_to_ECalDep, const int tot_layers, const int crit_layers);
  bool CurvatureResolutionFilter(int id, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, size_t>& ID_to_index, dunemc_plotting& plotting_vars, double pixel_spacing_cm);
  void EraseDescendants(int motherID, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID);
  bool IsParticleSelected(const int iSample, const int iEvent, const int iParticle);

  double pot;
  double _BeRPA_cvwgt = 1;

  //Geant vectors
  void clearBranchVectors();
  std::vector<float>* _MCVertX=0;
  std::vector<float>*  _MCVertY=0;
  std::vector<float>* _MCVertZ=0;
  std::vector<float>* _MCNuPx=0;
  std::vector<float>* _MCNuPy=0;
  std::vector<float>* _MCNuPz=0;
  std::vector<int>* _IsNC=0;
  std::vector<int>* _MCMode=0;
  std::vector<float> *_MCPStartX=0;
  std::vector<float> *_MCPStartY=0;
  std::vector<float> *_MCPStartZ=0;
  std::vector<float> *_MCPEndX=0;
  std::vector<float> *_MCPEndY=0;
  std::vector<float> *_MCPEndZ=0;
  std::vector<float> *_MCPStartPX=0;
  std::vector<float> *_MCPStartPY=0;
  std::vector<float> *_MCPStartPZ=0;
  std::vector<float> *_MCPEndPX=0;
  std::vector<float> *_MCPEndPY=0;
  std::vector<float> *_MCPEndPZ=0;
  std::vector<int> *_PDG = 0;
  std::vector<int> *_MCPTrkID=0;
  std::vector<std::string> *_MCPProc=0;
  std::vector<std::string> *_MCPEndProc=0;
  std::vector<int> *_MotherTrkID=0;
  std::vector<int> *_SimHitTrkID=0;
  std::vector<int> *_SimHitLayer=0;
  std::vector<float> *_SimHitEnergy=0;
  std::vector<float> *_SimHitX=0;
  std::vector<float> *_SimHitY=0;
  std::vector<float> *_SimHitZ=0;

  // Fixed detector variables
  double TPC_centre_x = -0.00012207;
  double TPC_centre_y = -150.473;
  double TPC_centre_z = 1486.;
  double BeamDirection[3] = {0.,-0.101,0.995};
  double X0 = 1193; //in cm From Federico's Kalman Filter Paper

  // MC-dependent detector variables
  double TPCInstrumentedLength;
  double TPCInstrumentedRadius;
  double ECALInnerRadius;
  double ECALOuterRadius;
  double ECALEndCapStart;
  double ECALEndCapEnd;

  // MC-independent configurable detector parameters
  double TPCFidLength;
  double TPCFidRadius;
  double B_field;
  double momentum_resolution_threshold;
  double pixel_spacing;
  double spatial_resolution;
  double adc_sampling_frequency;
  double drift_velocity;
  double downsampling;
  int crit_layers;
  double edepcrit_threshold;

  const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    {"TrueNeutrinoEnergy",kTrueNeutrinoEnergy},
    {"Mode",kMode},
    {"OscillationChannel",kOscChannel},
    {"TrueXPos",kTrueXPos},
    {"TrueYPos",kTrueYPos},
    {"TrueZPos",kTrueZPos},
    {"TrueRad",kTrueRad},
    {"TrueLepEnergy",kTrueLepEnergy},
    {"LepPT",kLepPT},
    {"LepPZ",kLepPZ},
    {"LepTheta",kLepTheta},
    {"LepPhi",kLepPhi},
    {"LepP",kLepP},
    {"LepTrackLengthYZ",kLepTrackLengthYZ},
    {"LepBAngle",kLepBAngle},
    {"TrueQ0",kTrueQ0},
    {"TrueQ3",kTrueQ3},
    {"Event_IsAccepted",kEvent_IsAccepted},
    {"InFDV",kInFDV},
    {"IsCC",kIsCC},
    {"EPi0",kEPi0},
    {"NPi0",kNPi0},
  };

  const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    {kTrueNeutrinoEnergy,"TrueNeutrinoEnergy"},
    {kMode,"Mode"},
    {kOscChannel,"OscillationChannel"},
    {kTrueXPos,"TrueXPos"},
    {kTrueYPos,"TrueYPos"},
    {kTrueZPos,"TrueZPos"},
    {kTrueRad,"TrueRad"},
    {kTrueLepEnergy,"TrueLepEnergy"},
    {kLepPT,"LepPT"},
    {kLepPZ,"LepPZ"},
    {kLepTheta,"LepTheta"},
    {kLepPhi,"LepPhi"},
    {kLepBAngle,"LepBAngle"},
    {kLepP,"LepP"},
    {kLepTrackLengthYZ,"LepTrackLengthYZ"},
    {kTrueQ0,"TrueQ0"},
    {kTrueQ3,"TrueQ3"},
    {kEvent_IsAccepted,"Event_IsAccepted"},
    {kInFDV,"InFDV"},
    {kIsCC,"IsCC"},
    {kEPi0,"EPi0"},
    {kNPi0,"NPi0"},
  };
    
  const std::unordered_map<std::string, int> KinematicVectorsDUNE = {
    {"Particle_EVis",kParticle_EVis},
    {"Particle_Momentum",kParticle_Momentum},
    {"Particle_EndMomentum",kParticle_EndMomentum},
    {"Particle_TransverseMomentum",kParticle_TransverseMomentum},
    {"Particle_BAngle",kParticle_BAngle},
    {"Particle_BeamAngle",kParticle_BeamAngle},
    {"Particle_IsAccepted",kParticle_IsAccepted},
    {"Particle_IsCurvatureResolved",kParticle_IsCurvatureResolved},
    {"Particle_IsDecayed",kParticle_IsDecayed},
    {"Particle_PDG",kParticle_PDG},
    {"Particle_IsStoppedInTPC",kParticle_IsStoppedInTPC},
    {"Particle_IsStoppedInECal",kParticle_IsStoppedInECal},
    {"Particle_IsStoppedInBarrel",kParticle_IsStoppedInBarrel},
    {"Particle_IsStoppedInEndCap",kParticle_IsStoppedInEndCap},
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
    {"Particle_EndDepth",kParticle_EndDepth},
    {"Particle_EndX",kParticle_EndX},
    {"Particle_EndY",kParticle_EndY},
    {"Particle_EndZ",kParticle_EndZ},
    {"Particle_StartX",kParticle_StartX},
    {"Particle_EDepCrit",kParticle_EDepCrit},
  };

  const std::unordered_map<int, std::string> ReversedKinematicVectorsDUNE = {
    {kParticle_EVis,"Particle_EVis"},
    {kParticle_Momentum,"Particle_Momentum"},
    {kParticle_EndMomentum,"Particle_EndMomentum"},
    {kParticle_TransverseMomentum,"Particle_TransverseMomentum"},
    {kParticle_BAngle,"Particle_BAngle"},
    {kParticle_BeamAngle,"Particle_BeamAngle"},
    {kParticle_IsAccepted,"Particle_IsAccepted"},
    {kParticle_IsCurvatureResolved,"Particle_IsCurvatureResolved"},
    {kParticle_IsDecayed,"Particle_IsDecayed"},
    {kParticle_PDG,"Particle_PDG"},
    {kParticle_IsStoppedInTPC,"Particle_IsStoppedInTPC"},
    {kParticle_IsStoppedInECal,"Particle_IsStoppedInECal"},
    {kParticle_IsStoppedInBarrel,"Particle_IsStoppedInBarrel"},
    {kParticle_IsStoppedInEndCap,"Particle_IsStoppedInEndCap"},
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
    {kParticle_EndDepth,"Particle_EndDepth"},
    {kParticle_EndX,"Particle_EndX"},
    {kParticle_EndY,"Particle_EndY"},
    {kParticle_EndZ,"Particle_EndZ"},
    {kParticle_StartX,"Particle_StartX"},
    {kParticle_EDepCrit,"Particle_EDepCrit"},
  };
    
};

#endif
