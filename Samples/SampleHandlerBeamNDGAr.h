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
  
  enum KinematicVecs {kPrim_EVis, kPrim_Momentum, kPrim_EndMomentum, kPrim_TransverseMomentum, 
    kPrim_BAngle, kPrim_BeamAngle, kPrim_IsAccepted, kPrim_IsCurvatureResolved, kPrim_IsDecayed, kPrim_PDG,
    kPrim_IsStoppedInTPC, kPrim_IsStoppedInECal, kPrim_IsStoppedInBarrel, kPrim_IsStoppedInEndCap, kPrim_IsStoppedInGap, 
    kPrim_IsStoppedInEndGap, kPrim_IsStoppedInBarrelGap, kPrim_IsEscaped, kPrim_NTurns, kPrim_NHits,
    kPrim_TrackLengthYZ, kPrim_MomResMS, kPrim_MomResYZ, kPrim_MomResX, kPrim_StartR2, kPrim_EndR, 
    kPrim_EndDepth, kPrim_EndX, kPrim_EndY, kPrim_EndZ, kPrim_StartX, kPrim_EDepCrit, kPrim_IsContained, kPrim_TPCEDepFrac,
    kShower_DCalBoundary, kShower_Energy, kShower_BAngle, kShower_IsContained, kShower_PDG, kShower_CosNorm,
    kPhoton_Energy, kPhoton_EndX, kPhoton_EndY, kPhoton_EndZ};

protected:
  //Functions required by core
  void Init() override;
  int SetupExperimentMC() override;
  void SetupFDMC() override;
  void AddAdditionalWeightPointers() override;
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
  std::vector<BeamNDSampleInfo> beamNDGArSampleDetails;

  //NDGAr-specific functions
  double FindNHits(double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature, double theta_start, double theta_spanned, int charge);
  bool isCoordOnTrack(int charge, double ycoord, double zcoord, double centre_circle_y, double centre_circle_z, double theta_start, double theta_spanned);
  double CalcBeta(double p_mag, double& bg, double& gamma, double pdgmass);
  int GetChargeFromPDG(int pdg);
  bool IsResolvedFromCurvature(dunemc_plotting& plotting_vars, size_t i_anapart, double pixel_spacing_cm);
  double GetCalDepth(double x, double y, double z);
  double GetDCalBoundary(const std::vector<double>& pos, const std::vector<double>& dir, size_t& boundary_index);
  double DepthToLayer(double depth, double r);
  double CalcEDepCal(int motherID, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, std::vector<double>>& ID_to_ECalDep, const int tot_layers);
  bool CurvatureResolutionFilter(int id, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, size_t>& ID_to_index, dunemc_plotting& plotting_vars, double pixel_spacing_cm);
  bool IsPrimContained(int id, const std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, size_t>& ID_to_index, 
                       const std::unordered_map<int, std::vector<double>>& eID_to_showerstart, const std::unordered_map<int, std::pair<double, double>>& pID_to_EDep,
                       dunemc_plotting& plotting_vars);
  void EraseDescendants(int motherID, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID);
  bool IsParticleSelected(const int iSample, const int iEvent, const int iParticle);
  void FillGeoVars();

  double _BeRPA_cvwgt = 1;
  
  // FastGArSim anatree inputs
  int _EventID;
  void clearBranchVectors();
  void fixCoordinates();
  std::vector<float> *_MCPStartX=nullptr;
  std::vector<float> *_MCPStartY=nullptr;
  std::vector<float> *_MCPStartZ=nullptr;
  std::vector<float> *_MCPEndX=nullptr;
  std::vector<float> *_MCPEndY=nullptr;
  std::vector<float> *_MCPEndZ=nullptr;
  std::vector<float> *_MCPStartPX=nullptr;
  std::vector<float> *_MCPStartPY=nullptr;
  std::vector<float> *_MCPStartPZ=nullptr;
  std::vector<float> *_MCPEndPX=nullptr;
  std::vector<float> *_MCPEndPY=nullptr;
  std::vector<float> *_MCPEndPZ=nullptr;
  std::vector<float> *_MCPCalPX=nullptr;
  std::vector<float> *_MCPCalPY=nullptr;
  std::vector<float> *_MCPCalPZ=nullptr;
  std::vector<int> *_MCPPDG=nullptr;
  std::vector<int> *_MCPTrkID=nullptr;
  std::vector<int> *_MCPMotherTrkID=nullptr;
  std::vector<int> *_TPCHitTrkID=nullptr;
  std::vector<float> *_TPCHitEnergy=nullptr;
  std::vector<float> *_TPCHitX=nullptr;
  std::vector<float> *_TPCHitY=nullptr;
  std::vector<float> *_TPCHitZ=nullptr;
  std::vector<bool> *_TPCHitIsSec=nullptr;
  std::vector<int> *_CalHitTrkID=nullptr;
  std::vector<int> *_CalHitLayer=nullptr;
  std::vector<float> *_CalHitEnergy=nullptr;
  std::vector<bool> *_CalHitIsSec=nullptr;
  std::vector<float> *_CalHitTime=nullptr;
  std::vector<float> *_CalHitX=nullptr;
  std::vector<float> *_CalHitY=nullptr;
  std::vector<float> *_CalHitZ=nullptr;
  std::vector<int> *_MuIDHitTrkID=nullptr;
  std::vector<float> *_MuIDHitEnergy=nullptr;
  std::vector<float> *_MuIDHitX=nullptr;
  std::vector<float> *_MuIDHitY=nullptr;
  std::vector<float> *_MuIDHitZ=nullptr;

  // FastGArSim geotree inputs
  double _TPCRad;
  double _TPCLen;
  double _BField;
  int _NumCalSides;
  double _BarrelGap;
  double _EndCapGap;
  double _HGAbsWidth;
  double _LGAbsWidth;
  double _HGSciWidth;
  double _LGSciWidth;
  double _HGBoardWidth;
  int _NBarrelHG;
  int _NBarrelLG;
  int _NEndCapHG;
  int _NEndCapLG;

  // Genie inputs
  double _Enu;
  double _PXnu;
  double _PYnu;
  double _PZnu;
  double _Elep;
  double _PXlep;
  double _PYlep;
  double _PZlep;
  int _nuPDG;
  bool _isCC;
  int _npip;
  int _npim;
  int _npi0;
  int _neut_code;

  // TPC dimensions
  double TPCFidLength;
  double TPCFidRadius;
  double TPCInstrumentedLength;
  double TPCInstrumentedRadius;
  double ECALInnerRadius;
  double ECALOuterRadius;
  double ECALEndCapStart;
  double ECALEndCapEnd;
  double ECALSciX0;
  std::vector<std::vector<double>> outerECalP;
  std::vector<std::vector<double>> outerECalA;
  
  double TPC_centre_x = 0.;
  double TPC_centre_y = 0.;
  double TPC_centre_z = 0.;
  double BeamDirection[3] = {0.,0.,1.};
  // double BeamDirection[3] = {0.,-0.101,0.995};

  double X0 = 1193; //in cm From Federico's Kalman Filter Paper

  //configurable detector parameters
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
    {"Prim_EVis",kPrim_EVis},
    {"Prim_Momentum",kPrim_Momentum},
    {"Prim_EndMomentum",kPrim_EndMomentum},
    {"Prim_TransverseMomentum",kPrim_TransverseMomentum},
    {"Prim_BAngle",kPrim_BAngle},
    {"Prim_BeamAngle",kPrim_BeamAngle},
    {"Prim_IsAccepted",kPrim_IsAccepted},
    {"Prim_IsCurvatureResolved",kPrim_IsCurvatureResolved},
    {"Prim_IsDecayed",kPrim_IsDecayed},
    {"Prim_PDG",kPrim_PDG},
    {"Prim_IsStoppedInTPC",kPrim_IsStoppedInTPC},
    {"Prim_IsStoppedInECal",kPrim_IsStoppedInECal},
    {"Prim_IsStoppedInBarrel",kPrim_IsStoppedInBarrel},
    {"Prim_IsStoppedInEndCap",kPrim_IsStoppedInEndCap},
    {"Prim_IsStoppedInGap",kPrim_IsStoppedInGap},
    {"Prim_IsStoppedInEndGap",kPrim_IsStoppedInEndGap},
    {"Prim_IsStoppedInBarrelGap",kPrim_IsStoppedInBarrelGap},
    {"Prim_IsEscaped",kPrim_IsEscaped},
    {"Prim_NTurns",kPrim_NTurns},
    {"Prim_NHits",kPrim_NHits},
    {"Prim_TrackLengthYZ",kPrim_TrackLengthYZ},
    {"Prim_MomResMS",kPrim_MomResMS},
    {"Prim_MomResYZ",kPrim_MomResYZ},
    {"Prim_MomResX",kPrim_MomResX},
    {"Prim_StartR2",kPrim_StartR2},
    {"Prim_EndR",kPrim_EndR},
    {"Prim_EndDepth",kPrim_EndDepth},
    {"Prim_EndX",kPrim_EndX},
    {"Prim_EndY",kPrim_EndY},
    {"Prim_EndZ",kPrim_EndZ},
    {"Prim_StartX",kPrim_StartX},
    {"Prim_EDepCrit",kPrim_EDepCrit},
    {"Prim_TPCEDepFrac",kPrim_TPCEDepFrac},
    {"Prim_IsContained",kPrim_IsContained},
    {"Shower_DCalBoundary",kShower_DCalBoundary},
    {"Shower_Energy",kShower_Energy},
    {"Shower_BAngle",kShower_BAngle},
    {"Shower_IsContained",kShower_IsContained},
    {"Shower_PDG",kShower_PDG},
    {"Shower_CosNorm",kShower_CosNorm},
    {"Photon_Energy",kPhoton_Energy},
    {"Photon_EndX",kPhoton_EndX},
    {"Photon_EndY",kPhoton_EndY},
    {"Photon_EndZ",kPhoton_EndZ},
  };

  const std::unordered_map<int, std::string> ReversedKinematicVectorsDUNE = {
    {kPrim_EVis,"Prim_EVis"},
    {kPrim_Momentum,"Prim_Momentum"},
    {kPrim_EndMomentum,"Prim_EndMomentum"},
    {kPrim_TransverseMomentum,"Prim_TransverseMomentum"},
    {kPrim_BAngle,"Prim_BAngle"},
    {kPrim_BeamAngle,"Prim_BeamAngle"},
    {kPrim_IsAccepted,"Prim_IsAccepted"},
    {kPrim_IsCurvatureResolved,"Prim_IsCurvatureResolved"},
    {kPrim_PDG,"Prim_PDG"},
    {kPrim_IsStoppedInTPC,"Prim_IsStoppedInTPC"},
    {kPrim_IsStoppedInECal,"Prim_IsStoppedInECal"},
    {kPrim_IsStoppedInBarrel,"Prim_IsStoppedInBarrel"},
    {kPrim_IsStoppedInEndCap,"Prim_IsStoppedInEndCap"},
    {kPrim_IsStoppedInGap,"Prim_IsStoppedInGap"},
    {kPrim_IsStoppedInEndGap,"Prim_IsStoppedInEndGap"},
    {kPrim_IsStoppedInBarrelGap,"Prim_IsStoppedInBarrelGap"},
    {kPrim_IsEscaped,"Prim_IsEscaped"},
    {kPrim_NTurns,"Prim_NTurns"},
    {kPrim_NHits,"Prim_NHits"},
    {kPrim_TrackLengthYZ,"Prim_TrackLengthYZ"},
    {kPrim_MomResMS,"Prim_MomResMS"},
    {kPrim_MomResYZ,"Prim_MomResYZ"},
    {kPrim_MomResX,"Prim_MomResX"},
    {kPrim_StartR2,"Prim_StartR2"},
    {kPrim_EndR,"Prim_EndR"},
    {kPrim_EndDepth,"Prim_EndDepth"},
    {kPrim_EndX,"Prim_EndX"},
    {kPrim_EndY,"Prim_EndY"},
    {kPrim_EndZ,"Prim_EndZ"},
    {kPrim_StartX,"Prim_StartX"},
    {kPrim_EDepCrit,"Prim_EDepCrit"},
    {kPrim_TPCEDepFrac,"Prim_TPCEDepFrac"},
    {kPrim_IsContained,"Prim_IsContained"},
    {kShower_DCalBoundary,"Shower_DCalBoundary"},
    {kShower_Energy,"Shower_Energy"},
    {kShower_BAngle,"Shower_BAngle"},
    {kShower_IsContained,"Shower_IsContained"},
    {kShower_PDG,"Shower_PDG"},
    {kShower_CosNorm,"Shower_CosNorm"},
    {kPhoton_Energy,"Photon_Energy"},
    {kPhoton_EndX,"Photon_EndX"},
    {kPhoton_EndY,"Photon_EndY"},
    {kPhoton_EndZ,"Photon_EndZ"},
  };
    
};

#endif
