#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"

#if defined(MaCh3_DUNE_USE_SRProxy) && (MaCh3_DUNE_USE_SRProxy==1)
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#endif

#pragma GCC diagnostic pop

SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerBase(mc_version_, xsec_cov_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  std::vector<std::string> EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);
  IsELike.resize(GetNSamples());
  
  ExposureScaling = Get<double>(SampleManager->raw()["AnalysisOptions"]["ExposureScaling"],__FILE__,__LINE__);
  for(int iSample=0;iSample<GetNSamples();iSample++){
    const std::string TempTitle = EnabledSamples[iSample];
    IsELike[iSample] = Get<int>(SampleManager->raw()[TempTitle]["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  }

  //DB Value used to determine selection criteria for FC and PC separation in function of 'walldist' variable
  FCPCSeparation = Get<double>(SampleManager->raw()["AnalysisOptions"]["FCPCSeparation"],__FILE__,__LINE__);

  //DB Define the names of the samples we're performing event selection for
  EventSelectionNames[kEventSel_FC_NuE]  = "FC_nueselec";
  EventSelectionNames[kEventSel_FC_NuMu] = "FC_numuselec";
  EventSelectionNames[kEventSel_FC_NC]   = "FC_ncselec";
  EventSelectionNames[kEventSel_PC_NuE]  = "PC_nueselec";
  EventSelectionNames[kEventSel_PC_NuMu] = "PC_numuselec";
  EventSelectionNames[kEventSel_PC_NC]   = "PC_ncselec";    

  //DB Create a map between these event selections to the defined samples in the config
  for (int iSelection=0;iSelection<nEventSelections;iSelection++) {
    EventSelection_to_SampleIndex_Map[iSelection] = kEventSel_Unknown;
    for (size_t iDefinedSample=0;iDefinedSample<SampleDetails.size();iDefinedSample++) {
      if (EventSelectionNames[iSelection] == SampleDetails[iDefinedSample].SampleTitle) {
	EventSelection_to_SampleIndex_Map[iSelection] = static_cast<int>(iDefinedSample); 
      }
    }
  }

  //DB Check atleast one of the event selections is in the config
  bool CheckVal = false;
  for (int iSelection=0;iSelection<nEventSelections;iSelection++) {
    if (EventSelection_to_SampleIndex_Map[iSelection] != kEventSel_Unknown) {
      CheckVal = true;
    }
  }
  if (CheckVal == false) {
    MACH3LOG_ERROR("No Event Selections match Defined Samples from Config");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
}

void SampleHandlerAtm::InititialiseData() {
  Reweight();
  for (int iSample = 0; iSample < GetNSamples(); iSample++) {
    AddData(iSample, GetMCArray(iSample));
  }
}

void SampleHandlerAtm::SetupSplines() {
  SplineHandler = nullptr;
}

void SampleHandlerAtm::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCEvents[i].total_weight_pointers.push_back(&(ExposureScaling));
  }  
}

int SampleHandlerAtm::SetupExperimentMC() {
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;  

  std::string InputFileName = Get<std::string>(SampleManager->raw()["InputFiles"]["FileName"],__FILE__,__LINE__);

  TFile *InputFile = TFile::Open(InputFileName.c_str(),"READ");
  if (!InputFile || InputFile->IsZombie()) {
    MACH3LOG_ERROR("Could not open input CAF file: {}",InputFileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TTree *cafTree, *weightsTree;
  InputFile->GetObject("cafTree",cafTree);
  if (!cafTree) {
    MACH3LOG_ERROR("Could not find cafTree in input CAF file: {}",InputFileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  InputFile->GetObject("weights",weightsTree);
  if (!weightsTree) {
    MACH3LOG_ERROR("Could not find weights tree in input CAF file: {}",InputFileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  double xsec_w, flux_nue_w, flux_numu_w;
  weightsTree->SetBranchAddress("xsec",&xsec_w);
  weightsTree->SetBranchAddress("flux_nue",&flux_nue_w);
  weightsTree->SetBranchAddress("flux_numu",&flux_numu_w);

#if defined(MaCh3_DUNE_USE_SRProxy) && (MaCh3_DUNE_USE_SRProxy==1)
  caf::StandardRecordProxy* sr = new caf::StandardRecordProxy(cafTree, "rec");    
#else  
  caf::StandardRecord* sr = new caf::StandardRecord();
  cafTree->SetBranchStatus("*", 1);
  cafTree->SetBranchAddress("rec", &sr);
#endif
  
  int nTreeEntries = static_cast<int>(cafTree->GetEntries());
  
  //================================================================================================

  for (int iTreeEntry=0;iTreeEntry<nTreeEntries;iTreeEntry++) {
    weightsTree->GetEntry(iTreeEntry);
    
#if defined(MaCh3_DUNE_USE_SRProxy) && (MaCh3_DUNE_USE_SRProxy==1)      
    cafTree->LoadTree(iTreeEntry);
#else
    cafTree->GetEntry(iTreeEntry);
#endif

    if ((iTreeEntry % (nTreeEntries/10))==0) {
      MACH3LOG_INFO("\tProcessing entry: {}/{}",iTreeEntry,nTreeEntries);
    }
    
    if(sr->common.ixn.pandora.size() != 1) {
      MACH3LOG_TRACE("Skipping entry {}/{} -> Number of neutrino slices found in event: {}",iTreeEntry,nTreeEntries,sr->common.ixn.pandora.size());
      continue;
    }

    /*
    int RunNumber = sr->meta.fd_hd.run;
    int SubRunNumber = sr->meta.fd_hd.subrun;
    int EventNumber = sr->meta.fd_hd.event;
    */
    
    std::vector<double> CVNScores = std::vector<double>(nCVN_Scores);
    CVNScores[kCVN_NuE] = sr->common.ixn.pandora[0].nuhyp.cvn.nue;
    CVNScores[kCVN_NuMu] = sr->common.ixn.pandora[0].nuhyp.cvn.numu;
    CVNScores[kCVN_NC] = sr->common.ixn.pandora[0].nuhyp.cvn.nc;    

    //DB Theoretical max value should be ShortestDimensioninAV/2 which is less than 1e4 -- so dummy value of 1e8 is fine
    double MinDistToWall = 1e8;
    for (size_t iPart=0;iPart<sr->common.ixn.pandora[0].part.pandora.size();iPart++) {
      //DB Info from PG -- ignore any hits associated with HitCollection classified objects
      if (sr->common.ixn.pandora[0].part.pandora[iPart].origRecoObjType == caf::RecoObjType::kHitCollection) {continue;}
      if (sr->common.ixn.pandora[0].part.pandora[iPart].walldist < MinDistToWall) {MinDistToWall = sr->common.ixn.pandora[0].part.pandora[iPart].walldist;}
    }
    
    int SampleIndex = ReturnSampleIdentifier(CVNScores, MinDistToWall);
    if (SampleIndex == kEventSel_Unknown) {
      continue;
    }

    TVector3 RecoNuMomentumVector;
    double RecoENu;
    if (IsELike[SampleIndex]) {
      RecoENu = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.x,sr->common.ixn.pandora[0].dir.heshw.y,sr->common.ixn.pandora[0].dir.heshw.z)).Unit();
    } else {
      RecoENu = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.x,sr->common.ixn.pandora[0].dir.lngtrk.y,sr->common.ixn.pandora[0].dir.lngtrk.z)).Unit();      
    }
    double RecoCZ = -RecoNuMomentumVector.y(); // +Y in CAF files translates to +Z in typical CosZ
    if (std::isnan(RecoCZ)) {
      MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Cosine Z is NAN",iTreeEntry,nTreeEntries);
      continue;
    }
    if (std::isnan(RecoENu)) {
      MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Neutrino Energy is NAN",iTreeEntry,nTreeEntries);
      continue;
    }

    auto& OscillationChannels = SampleDetails[SampleIndex].OscChannels;    
    int InteractingPDG = sr->mc.nu[0].pdg;

    int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
    if (!sr->mc.nu[0].iscc) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon

    double TrueNeutrinoEnergy = static_cast<double>(sr->mc.nu[0].E);
    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.x,sr->mc.nu[0].momentum.y,sr->mc.nu[0].momentum.z)).Unit();

    struct dunemc_atm currentEvent_FromNuE;
    
    currentEvent_FromNuE.rw_erec = RecoENu;
    currentEvent_FromNuE.rw_theta = RecoCZ;
    currentEvent_FromNuE.SampleIndex = SampleIndex;
    currentEvent_FromNuE.nupdg = InteractingPDG;
    currentEvent_FromNuE.nupdgUnosc = (InteractingPDG > 0) ? 12 : -12;
    currentEvent_FromNuE.OscChannelIndex = static_cast<double>(GetOscChannel(OscillationChannels, currentEvent_FromNuE.nupdgUnosc, currentEvent_FromNuE.nupdg));
    currentEvent_FromNuE.mode = M3Mode;
    currentEvent_FromNuE.rw_isCC = sr->mc.nu[0].iscc;
    currentEvent_FromNuE.Target = kTarget_Ar;
    currentEvent_FromNuE.enu_true = TrueNeutrinoEnergy;
    currentEvent_FromNuE.coszenith_true = -TrueNuMomentumVector.y(); // +Y in CAF files translates to +Z in typical CosZ
    currentEvent_FromNuE.flux_w = xsec_w*flux_nue_w;
    currentEvent_FromNuE.MinDistToWall = MinDistToWall;
    
    struct dunemc_atm currentEvent_FromNuMu = currentEvent_FromNuE;
    
    currentEvent_FromNuMu.nupdgUnosc = (InteractingPDG > 0) ? 14 : -14;
    currentEvent_FromNuMu.OscChannelIndex = static_cast<double>(GetOscChannel(OscillationChannels, currentEvent_FromNuMu.nupdgUnosc, currentEvent_FromNuMu.nupdg));
    currentEvent_FromNuMu.flux_w = xsec_w*flux_numu_w;

    dunemcSamples.emplace_back(std::move(currentEvent_FromNuE));
    dunemcSamples.emplace_back(std::move(currentEvent_FromNuMu));    
  }

  //================================================================================================
  gErrorIgnoreLevel = CurrErrorLevel;

  delete sr;
  delete cafTree;
  delete weightsTree;
  delete InputFile;
  
  return static_cast<int>(dunemcSamples.size());
}

void SampleHandlerAtm::SetupMC() {
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent) {
    MCEvents[iEvent].enu_true = dunemcSamples[iEvent].enu_true;
    MCEvents[iEvent].isNC = !dunemcSamples[iEvent].rw_isCC;
    MCEvents[iEvent].nupdg = dunemcSamples[iEvent].nupdg;
    MCEvents[iEvent].nupdgUnosc = dunemcSamples[iEvent].nupdgUnosc;
    MCEvents[iEvent].NominalSample = dunemcSamples[iEvent].SampleIndex;
    
    MCEvents[iEvent].coszenith_true = dunemcSamples[iEvent].coszenith_true;
  }
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(const int KinPar, int iEvent) const {
  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    return &(dunemcSamples[iEvent].enu_true);
  case kRecoNeutrinoEnergy:
    return &(dunemcSamples[iEvent].rw_erec);
  case kTrueCosZ:
    return &(dunemcSamples[iEvent].coszenith_true);
  case kRecoCosZ:
    return &(dunemcSamples[iEvent].rw_theta);
  case kOscChannel:
    return &(dunemcSamples[iEvent].OscChannelIndex);
  case kMode:
    return &(dunemcSamples[iEvent].mode);
  case kTargetNucleus:
    return &(dunemcSamples[iEvent].Target);
  case kMinDistToWall:
    return &(dunemcSamples[iEvent].MinDistToWall);
  default:
    MACH3LOG_ERROR("Unknown KinPar: {}",static_cast<int>(KinPar));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}


double SampleHandlerAtm::ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

int SampleHandlerAtm::ReturnSampleIdentifier(std::vector<double> CVNScores, double MinDistanceToWall) {
  bool IsFullyContained = false;
  
  if (MinDistanceToWall > 1e4 || MinDistanceToWall < 0) { //DB: ToDo Work out theoretical maximum
    return kEventSel_Unknown;
  } else if (MinDistanceToWall > FCPCSeparation) {
    IsFullyContained = true;
  } else {
    IsFullyContained = false;
  }

  //int EventSelection = static_cast<int>(std::distance(CVNScores.begin(), max_element(CVNScores.begin(), CVNScores.end())));
  //DB Not making the 0.55 and 0.56 magic numbers config-read because will eventually move to the argmax style selection (commented out on line above...)
  int EventSelection = kCVN_NC;
  if (CVNScores[kCVN_NuMu] > 0.56) { 
    EventSelection = kCVN_NuMu;
  } else if (CVNScores[kCVN_NuE] > 0.55) {
    EventSelection = kCVN_NuE;
  }

  int SampleIndex = kEventSel_Unknown;
  if (IsFullyContained) {
    if (EventSelection == kCVN_NuE)  {SampleIndex = EventSelection_to_SampleIndex_Map[kEventSel_FC_NuE]; }
    if (EventSelection == kCVN_NuMu) {SampleIndex = EventSelection_to_SampleIndex_Map[kEventSel_FC_NuMu];}
    if (EventSelection == kCVN_NC)   {SampleIndex = EventSelection_to_SampleIndex_Map[kEventSel_FC_NC];  }    
  } else {
    if (EventSelection == kCVN_NuE)  {SampleIndex = EventSelection_to_SampleIndex_Map[kEventSel_PC_NuE]; }
    if (EventSelection == kCVN_NuMu) {SampleIndex = EventSelection_to_SampleIndex_Map[kEventSel_PC_NuMu];}
    if (EventSelection == kCVN_NC)   {SampleIndex = EventSelection_to_SampleIndex_Map[kEventSel_PC_NC];  }    
  }
  
  return SampleIndex;
}
