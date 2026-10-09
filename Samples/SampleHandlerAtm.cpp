#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"

//Used to read from CAF files
#if defined(MaCh3_DUNE_USE_SRProxy) && (MaCh3_DUNE_USE_SRProxy==1)
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#endif
#include "duneanaobj/StandardRecord/StandardRecord.h"

//Used to loop through ROOT files in a directory
#include <filesystem>

//Used to read from an Eigen binary blob
#include <Eigen/Dense>
#include "DUNEUtils.h"

#pragma GCC diagnostic pop

SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerBase(mc_version_, xsec_cov_, Oscillator_) {
  //DB Register the kinematic parameter maps
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  //DB IsELike determines which reconstruction algorithm gets applied to which samples
  std::vector<std::string> EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);
  IsELike.resize(GetNSamples());
  for(int iSample=0;iSample<GetNSamples();iSample++){
    const std::string TempTitle = EnabledSamples[iSample];
    IsELike[iSample] = Get<int>(SampleManager->raw()[TempTitle]["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  }

  //DB Multiplicative scaling factor on the base exposure of the MC
  ExposureScaling = Get<double>(SampleManager->raw()["AnalysisOptions"]["ExposureScaling"],__FILE__,__LINE__);
 
  //DB Value used to determine selection criteria for FC and PC separation in function of 'walldist' variable
  FCPCSeparation = Get<double>(SampleManager->raw()["AnalysisOptions"]["FCPCSeparation"],__FILE__,__LINE__);

  //=============================================================
  //DB MC read in information
  //   Highest priority: If EigenInputFile is provided, MC is read in from this Eigen binary
  //   Middle priority: If InputFileName is provided, MC will be read in from this single CAF file
  //   Lowest priority: All CAF ROOT files in the 'InputFileDirectory' directory will be read in -- be warned
  
  if (SampleManager->raw()["InputFiles"]["EigenFile"]) {
    EigenInputFile = Get<std::string>(SampleManager->raw()["InputFiles"]["EigenFile"],__FILE__,__LINE__);
    EigenInputFileMD5Sum = Get<std::string>(SampleManager->raw()["InputFiles"]["EigenFileMD5Sum"],__FILE__,__LINE__);
  } else {
    EigenInputFile = "";
  }

  InputFileName = Get<std::string>(SampleManager->raw()["InputFiles"]["FileName"],__FILE__,__LINE__);
  InputFileDirectory = Get<std::string>(SampleManager->raw()["InputFiles"]["FileDirectory"],__FILE__,__LINE__);
  //=============================================================
  
  // Per-event spline weights (optional)
  if (SampleManager->raw()["InputFiles"]["InputSplines"]) {
    fInputSplines = Get<std::string>(SampleManager->raw()["InputFiles"]["InputSplines"],__FILE__,__LINE__);
  } else {
    fInputSplines = "";
  }

  //DB Input file containing the detector systematic spectra ratios, with the expectation that HistogramTitle is formatted as SYSTEMATICNAME_SAMPLENAME
  if (SampleManager->raw()["InputFiles"]["DetectorSystematicsFileName"]) {
    fDetectorSystematicsFileName = Get<std::string>(SampleManager->raw()["InputFiles"]["DetectorSystematicsFileName"],__FILE__,__LINE__);
  } else {
    fDetectorSystematicsFileName = "";
  }
  
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
  //DB By default, set the data to be the "nominal" MC spectrum - This can be overwritten at executable level
  Reweight();
  for (int iSample = 0; iSample < GetNSamples(); iSample++) {
    AddData(iSample, GetMCArray(iSample));
  }
}

void SampleHandlerAtm::SetupSplines() {
  //DB Because the event-by-event spline weight load-in is not "easily" supported by Core, do manual load-in using DUNE tools
  if(ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline) > 0) {
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    auto SplineFactory = SplineHandlerFactoryDUNE(ParHandler, Modes.get(), dunemcSamples, fInputSplines, SampleHandlerName);

    SplineHandler = std::move(SplineFactory.GetSplineHandler());
    if (SplineFactory.GetSplineType() == kBinned) {
      InitialiseSplineObject(); //Running the "normal" initialisation for binned splines
    } else if (SplineFactory.GetSplineType() == kMonolith) {
      InitialiseSplineObjectPerEvent(); //Running the per-event initialisation
    } else {
      MACH3LOG_ERROR("Unknown spline type found when setting up splines for sample {}", SampleHandlerName);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
  } else{
    MACH3LOG_INFO("Found no spline for this sample so I will not load or evaluate splines");
    SplineHandler = nullptr;
  }
}

void SampleHandlerAtm::InitialiseSplineObjectPerEvent(){
  //DB Do the manual hand-holding of intialising the event-by-event spline and setting the event weights
  auto SplineHandlerDUNE = dynamic_cast<MonolithSplineHandlerDUNE*>(SplineHandler.get());
  
  if (!SplineHandlerDUNE) {
    MACH3LOG_ERROR("SplineHandler is not of type MonolithSplineHandlerDUNE for sample {}. Cannot call InitialiseSplineObjectPerEvent.", SampleHandlerName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(SplineHandlerDUNE->RetPointer(static_cast<int>(i)));
  }
}

void SampleHandlerAtm::RegisterFunctionalParameters() {
  if (!ParHandler) {return;}

  //DB Throw all the detector systematics into a single weight -- calculate the detector systematics weight as (1.0 + ParVal*SpectraRatio)
  if (DetectorSystematicParameterNames.size()!=0) {
    RegisterIndividualFunctionalParameter(dunemcSamples, DetectorSystematicParameterNames, [](std::vector<double> const &ParVals, dunemc_atm &Event) {
      Event.TotalDetectorSystematicWeight = 1.0;
      
      for (size_t iSyst=0;iSyst<ParVals.size();iSyst++) {
	//DB The minus accounts for the fact we only want the "difference" in ratio
	Event.TotalDetectorSystematicWeight *= 1.0 + ParVals[iSyst]*(Event.DetectorSystematicRatios[iSyst] - 1.0);
      }
    });
  }
}

void SampleHandlerAtm::AddAdditionalWeightPointers() {
  //DB ToDo: Try to "group" as many fixed weights as possible to reduce needless multiple multiplications
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCEvents[i].total_weight_pointers.push_back(&(ExposureScaling));
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].TotalDetectorSystematicWeight));    
  }  
}

int SampleHandlerAtm::ReadFromEigen() {
  //DB Read the MC info from the Eigen binary blob -- but first check that the MD5 Checksum matches the expectation
  
  Eigen::MatrixXd Matrix;

  MACH3LOG_INFO("Reading Eigen matrix from:{}",EigenInputFile);
  CompareMD5Sum(EigenInputFileMD5Sum,EigenInputFile);
  ReadEigenMatrixFromFile(EigenInputFile.c_str(),Matrix);

  int nEntries = static_cast<int>(Matrix.rows());
  dunemcSamples.resize(nEntries);

  //DB The below needs modifying any time a new variable is read from, or derived from, the CAF file
  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    dunemcSamples[iEvent].rw_erec = Matrix(iEvent,Variables::rw_erec);
    dunemcSamples[iEvent].rw_ehad = Matrix(iEvent,Variables::rw_ehad);
    dunemcSamples[iEvent].rw_elep = Matrix(iEvent,Variables::rw_elep);
    dunemcSamples[iEvent].rw_theta = Matrix(iEvent,Variables::rw_theta);    
    dunemcSamples[iEvent].SampleIndex = static_cast<int>(Matrix(iEvent,Variables::SampleIndex));
    dunemcSamples[iEvent].nupdg = static_cast<int>(Matrix(iEvent,Variables::nupdg));
    dunemcSamples[iEvent].nupdgUnosc = static_cast<int>(Matrix(iEvent,Variables::nupdgUnosc));
    dunemcSamples[iEvent].OscChannelIndex = Matrix(iEvent,Variables::OscChannelIndex);
    dunemcSamples[iEvent].mode = Matrix(iEvent,Variables::mode);
    dunemcSamples[iEvent].rw_isCC = static_cast<int>(Matrix(iEvent,Variables::rw_isCC));
    dunemcSamples[iEvent].Target = static_cast<int>(Matrix(iEvent,Variables::Target));
    dunemcSamples[iEvent].enu_true = Matrix(iEvent,Variables::enu_true);
    dunemcSamples[iEvent].coszenith_true = Matrix(iEvent,Variables::coszenith_true);    
    dunemcSamples[iEvent].flux_w = Matrix(iEvent,Variables::flux_w);
    dunemcSamples[iEvent].MinDistToWall = Matrix(iEvent,Variables::MinDistToWall);
    dunemcSamples[iEvent].eid = static_cast<uint>(Matrix(iEvent,Variables::eid));
  }

  return nEntries;
}

void SampleHandlerAtm::TransferToEigen(std::string FileName) {
  Eigen::MatrixXd Matrix = Eigen::MatrixXd(dunemcSamples.size(),Variables::nVariables);

  //DB The below needs modifying any time a new variable is read from, or derived from, the CAF file
  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    Matrix(iEvent,Variables::rw_erec) = dunemcSamples[iEvent].rw_erec;
    Matrix(iEvent,Variables::rw_ehad) = dunemcSamples[iEvent].rw_ehad;
    Matrix(iEvent,Variables::rw_elep) = dunemcSamples[iEvent].rw_elep;
    Matrix(iEvent,Variables::rw_theta) = dunemcSamples[iEvent].rw_theta;
    Matrix(iEvent,Variables::SampleIndex) = dunemcSamples[iEvent].SampleIndex;
    Matrix(iEvent,Variables::nupdg) = dunemcSamples[iEvent].nupdg;
    Matrix(iEvent,Variables::nupdgUnosc) = dunemcSamples[iEvent].nupdgUnosc;
    Matrix(iEvent,Variables::OscChannelIndex) = dunemcSamples[iEvent].OscChannelIndex;
    Matrix(iEvent,Variables::mode) = dunemcSamples[iEvent].mode;
    Matrix(iEvent,Variables::rw_isCC) = dunemcSamples[iEvent].rw_isCC;    
    Matrix(iEvent,Variables::Target) = dunemcSamples[iEvent].Target;
    Matrix(iEvent,Variables::enu_true) = dunemcSamples[iEvent].enu_true;
    Matrix(iEvent,Variables::coszenith_true) = dunemcSamples[iEvent].coszenith_true;    
    Matrix(iEvent,Variables::flux_w) = dunemcSamples[iEvent].flux_w;
    Matrix(iEvent,Variables::MinDistToWall) = dunemcSamples[iEvent].MinDistToWall;
    Matrix(iEvent,Variables::eid) = dunemcSamples[iEvent].eid;        
  }

  MACH3LOG_INFO("Writing Eigen matrix to:{}",FileName);
  WriteEigenMatrixToFile(FileName.c_str(),Matrix);
}

int SampleHandlerAtm::SetupExperimentMC() {
  // DB If the Eigen input file is define, read from that. Otherwise read from CAF files
  if (EigenInputFile != "") {
    ReadFromEigen();
  } else {
    
    int CurrErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;  

    //DB The TChain object which contains the CAF file information
    TChain* cafTree = new TChain("cafTree");
    //DB The supplemental file that contains the post-simulation processed weights (e.g. flux)
    TChain* weightsTree = new TChain("weights");

    if (InputFileName != "") {
    //If reading a single CAF file      
      cafTree->Add(InputFileName.c_str());
      weightsTree->Add(InputFileName.c_str());
    } else if (InputFileDirectory != "") {
      //If reading multiple ROOT CAF files
      std::filesystem::path directory = InputFileDirectory;
      
      for (const auto& entry : std::filesystem::directory_iterator(directory)) {
	if (entry.is_regular_file() && entry.path().extension()==".root") {
	  std::string FileName = entry.path();
	  MACH3LOG_INFO("Adding file:{}",FileName);
	  cafTree->Add(FileName.c_str());
	  weightsTree->Add(FileName.c_str());
	}
      }
    } else {
      MACH3LOG_ERROR("Bad configuration - Neither Eigen binary, single CAF file, or file directory passed");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    //DB Grab the weight information from the post-processed weight information (e.g. flux)
    double xsec_w, flux_nue_w, flux_numu_w;
    weightsTree->SetBranchAddress("xsec",&xsec_w);
    weightsTree->SetBranchAddress("flux_nue",&flux_nue_w);
    weightsTree->SetBranchAddress("flux_numu",&flux_numu_w);
    
#if defined(MaCh3_DUNE_USE_SRProxy) && (MaCh3_DUNE_USE_SRProxy==1)
    MACH3LOG_INFO("Using the Standard Record Proxy reader");
    caf::StandardRecordProxy* sr = new caf::StandardRecordProxy(cafTree, "rec");

    //DB Need to know the offsets to deal with the manual file changing within the chain
    int currentTreeNumber = 0;
    Long64_t *treeOffsets = cafTree->GetTreeOffset();
    int nbTrees = cafTree->GetTreeOffsetLen();
#else
    MACH3LOG_INFO("Using the usual Standard Record reader");
    caf::StandardRecord* sr = new caf::StandardRecord();
    cafTree->SetBranchStatus("*", 1);
    cafTree->SetBranchAddress("rec", &sr);
#endif
  
    int nTreeEntries = static_cast<int>(cafTree->GetEntries());
    dunemcSamples.reserve(2 * nTreeEntries);    
    //================================================================================================
    //DB Read in the CAF file information
    
    for (int iTreeEntry=0;iTreeEntry<nTreeEntries;iTreeEntry++) {
      weightsTree->GetEntry(iTreeEntry);
      
#if defined(MaCh3_DUNE_USE_SRProxy) && (MaCh3_DUNE_USE_SRProxy==1)
      cafTree->LoadTree(iTreeEntry);

      if (currentTreeNumber < nbTrees - 1 && iTreeEntry == treeOffsets[currentTreeNumber+1]) {
	//DB We are changing tree and due to the inability of SRProxy to handle it correctly, we do it manually
	currentTreeNumber++;
	delete sr;
	sr = new caf::StandardRecordProxy(cafTree->GetTree(), "rec");
      }
#else
      cafTree->GetEntry(iTreeEntry);
#endif
      
      if ((iTreeEntry % (nTreeEntries/20))==0) {
	MACH3LOG_INFO("\tProcessing entry: {}/{}",iTreeEntry,nTreeEntries);
      }

      //DB Need exactly one reconstructed object for an event
      if(sr->common.ixn.pandora.size() != 1) {
	MACH3LOG_TRACE("Skipping entry {}/{} -> Number of neutrino slices found in event: {}",iTreeEntry,nTreeEntries,sr->common.ixn.pandora.size());
	continue;
      }
      
      /*
      //DB MC identification information
      int RunNumber = sr->meta.fd_hd.run;
      int SubRunNumber = sr->meta.fd_hd.subrun;
      int EventNumber = sr->meta.fd_hd.event;
      */
      
      std::vector<double> CVNScores = std::vector<double>(nCVN_Scores);
      CVNScores[kCVN_NuE] = sr->common.ixn.pandora[0].nuhyp.cvn.nue;
      CVNScores[kCVN_NuMu] = sr->common.ixn.pandora[0].nuhyp.cvn.numu;
      CVNScores[kCVN_NC] = sr->common.ixn.pandora[0].nuhyp.cvn.nc;    
      
      //PG Pre-check: If this event cannot pass selection cuts under either Fully Contained or Partially Contained assumption, skip it!
      //DB ToDo: Consider how to make this future-proof/less error prone
      if (ReturnSampleIdentifier(CVNScores, 1e8) == kEventSel_Unknown &&
	  ReturnSampleIdentifier(CVNScores, 0.0) == kEventSel_Unknown) {
	continue;
      }
      
      //DB Now, since the event has a valid candidate classification, we lazy-load and calculate MinDist:
      double MinDist = 1e8;
      auto const& PandoraParticles = sr->common.ixn.pandora[0].part.pandora;
      size_t nParts = PandoraParticles.size();
      for (size_t iPart=0;iPart<nParts;iPart++) {
	auto const& Part = PandoraParticles[iPart];
	//DB Info from PG -- ignore any hits associated with HitCollection classified objects
	if (Part.origRecoObjType == caf::RecoObjType::kHitCollection) {continue;}
	if (Part.walldist < MinDist) {MinDist = Part.walldist;}
      }

      //DB Check which sample the event belongs to based on the information read in
      int SampIndex = ReturnSampleIdentifier(CVNScores, MinDist);
      if (SampIndex == kEventSel_Unknown) {
	continue;
      }
      
      TVector3 RecoNuMomentumVector;
      double RecoENu;
      double RecoEHad;
      double RecoELep;

      //DB Calculate the reconstructed information depending on the requested reconstructed algorithms
      if (IsELike[SampIndex]) {
	RecoENu = sr->common.ixn.pandora[0].Enu.e_calo;
	RecoEHad = sr->common.ixn.pandora[0].Enu.e_had;
	RecoELep = RecoENu-RecoEHad;
	RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.x,sr->common.ixn.pandora[0].dir.heshw.y,sr->common.ixn.pandora[0].dir.heshw.z)).Unit();
      } else {
	RecoENu = sr->common.ixn.pandora[0].Enu.lep_calo;
	RecoEHad = sr->common.ixn.pandora[0].Enu.e_had;
	RecoELep = RecoENu-RecoEHad;	
	RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.x,sr->common.ixn.pandora[0].dir.lngtrk.y,sr->common.ixn.pandora[0].dir.lngtrk.z)).Unit();      
      }
      double RecoCZ = -RecoNuMomentumVector.y(); // +Y in CAF files translates to +Z in typical CosZ

      //DB Check if all the values are sensible
      if (std::isnan(RecoCZ)) {
	MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Cosine Z is NAN",iTreeEntry,nTreeEntries);
	continue;
      }
      if (std::isnan(RecoENu)) {
	MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Neutrino Energy is NAN",iTreeEntry,nTreeEntries);
	continue;
      }
      if (std::isnan(RecoEHad)) {
	MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Hadronic Energy is NAN",iTreeEntry,nTreeEntries);
	continue;
      }
      if (std::isnan(RecoELep)) {
	MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Leptonic Energy is NAN",iTreeEntry,nTreeEntries);
	continue;
      }            

      //DB Calculate the true kinematic information
      double TrueNeutrinoEnergy = static_cast<double>(sr->mc.nu[0].E);
      TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.x,sr->mc.nu[0].momentum.y,sr->mc.nu[0].momentum.z)).Unit();
      
      //DB Determine the PDG, oscillation channel, and MaCh3 index for the events GENIE mode
      auto& OscillationChannels = SampleDetails[SampIndex].OscChannels;    
      int InteractingPDG = sr->mc.nu[0].pdg;
      
      int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
      if (!sr->mc.nu[0].iscc) M3Mode += 14; //Account for no ability to distinguish CC/NC
      if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
      
      //DB Save the event information twice; once under the assumption of generated nue(bar) and the other generated as a numu(bar)
      struct dunemc_atm currentEvent_FromNuE;
      
      currentEvent_FromNuE.rw_erec = RecoENu;
      currentEvent_FromNuE.rw_ehad = RecoEHad;
      currentEvent_FromNuE.rw_elep = RecoELep;            
      currentEvent_FromNuE.rw_theta = RecoCZ;
      currentEvent_FromNuE.SampleIndex = SampIndex;
      currentEvent_FromNuE.nupdg = InteractingPDG;
      currentEvent_FromNuE.nupdgUnosc = (InteractingPDG > 0) ? 12 : -12;
      currentEvent_FromNuE.OscChannelIndex = static_cast<double>(GetOscChannel(OscillationChannels, currentEvent_FromNuE.nupdgUnosc, currentEvent_FromNuE.nupdg));
      currentEvent_FromNuE.mode = M3Mode;
      currentEvent_FromNuE.rw_isCC = sr->mc.nu[0].iscc;
      currentEvent_FromNuE.Target = kTarget_Ar;
      currentEvent_FromNuE.enu_true = TrueNeutrinoEnergy;
      currentEvent_FromNuE.coszenith_true = -TrueNuMomentumVector.y(); // +Y in CAF files translates to +Z in typical CosZ
      currentEvent_FromNuE.flux_w = xsec_w*flux_nue_w;
      currentEvent_FromNuE.MinDistToWall = MinDist;
      currentEvent_FromNuE.eid = static_cast<uint>(iTreeEntry);
      
      struct dunemc_atm currentEvent_FromNuMu = currentEvent_FromNuE;
      
      currentEvent_FromNuMu.nupdgUnosc = (InteractingPDG > 0) ? 14 : -14;
      currentEvent_FromNuMu.OscChannelIndex = static_cast<double>(GetOscChannel(OscillationChannels, currentEvent_FromNuMu.nupdgUnosc, currentEvent_FromNuMu.nupdg));
      currentEvent_FromNuMu.flux_w = xsec_w*flux_numu_w;

      //DB Save the information into the MC struct
      dunemcSamples.emplace_back(std::move(currentEvent_FromNuE));
      dunemcSamples.emplace_back(std::move(currentEvent_FromNuMu));    
    }
    
    //================================================================================================
    gErrorIgnoreLevel = CurrErrorLevel;
    
#if defined(MaCh3_DUNE_USE_SRProxy) && (MaCh3_DUNE_USE_SRProxy==1)  
    //PG Need to clear that static vector to avoid double free errors when exiting the program
    caf::SRBranchRegistry::clear();
#endif
    
    delete sr;
    delete cafTree;
    delete weightsTree;
  }

  //DB Regardless of however the information is read in, load the detector systematic ratio factors
  SetupDetectorSystematicRatios();

  //DB Return the number of events loaded
  return static_cast<int>(dunemcSamples.size());
}

void SampleHandlerAtm::SetupMC() {
  //DB Pass the Core-required information over 
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent) {
    MCEvents[iEvent].enu_true = dunemcSamples[iEvent].enu_true;
    MCEvents[iEvent].isNC = !dunemcSamples[iEvent].rw_isCC;
    MCEvents[iEvent].nupdg = dunemcSamples[iEvent].nupdg;
    MCEvents[iEvent].nupdgUnosc = dunemcSamples[iEvent].nupdgUnosc;
    MCEvents[iEvent].NominalSample = dunemcSamples[iEvent].SampleIndex;
    MCEvents[iEvent].coszenith_true = dunemcSamples[iEvent].coszenith_true;
  }
}

void SampleHandlerAtm::SetupDetectorSystematicRatios() {
  //DB Firstly set detector systematic weight to 1.0
  for (size_t iEvent=0;iEvent<dunemcSamples.size();iEvent++) {
    dunemcSamples[iEvent].TotalDetectorSystematicWeight = 1.0;
  }
  
  //DB Figure out whether we have detector systematics configured -- required that any Detector Systematic parameter's group is "DetectorSystematic"
  int NParams = ParHandler->GetNParameters();
  for (int iParam=0;iParam<NParams;iParam++) {
    if (ParHandler->IsParFromGroup(iParam,"DetectorSystematic")) {
      std::string ParamName = ParHandler->GetParFancyName(iParam);
      DetectorSystematicParameterNames.push_back(ParamName);
    }
  }

  //DB Check that both the Detector Systematic input file has been specified and that the parameters have been configured
  if (DetectorSystematicParameterNames.size()==0 && fDetectorSystematicsFileName!="") {
    MACH3LOG_ERROR("Detector Systematic Input file provided but no detector systematics (Group=\"DetectorSystematic\") found");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (DetectorSystematicParameterNames.size()!=0 && fDetectorSystematicsFileName=="") {
    MACH3LOG_ERROR("Detector systematics configured (Group=\"DetectorSystematic\") but no input file provided");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  //DB If no parameters are configured, we can return early
  if (DetectorSystematicParameterNames.size()==0) {
    return;
  }

  //DB Lets grab the Detector Systematic input file
  TFile* DetectorSystematicsFile = new TFile(fDetectorSystematicsFileName.c_str());
  if (!DetectorSystematicsFile || DetectorSystematicsFile->IsZombie()) {
    MACH3LOG_ERROR("Could not find file:{}",fDetectorSystematicsFileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //DB Load histograms with title SystematicName_SampleName into vector indexed [Systematic][Sample]
  std::vector<std::vector<TH1*>> RatioHistograms;
  RatioHistograms.resize(DetectorSystematicParameterNames.size());
  for (size_t iRatioHist=0;iRatioHist<RatioHistograms.size();iRatioHist++) {
    RatioHistograms[iRatioHist].resize(SampleDetails.size());
  }

  //DB Push the histograms into the vector
  for (size_t iParam=0;iParam<DetectorSystematicParameterNames.size();iParam++) {
    for (size_t iSample=0;iSample<SampleDetails.size();iSample++) {
      std::string ExpectedHistogramName = DetectorSystematicParameterNames[iParam]+"_"+SampleDetails[iSample].SampleTitle;
      TH1* Histogram = DetectorSystematicsFile->Get<TH1>(ExpectedHistogramName.c_str());
      if (!Histogram) {
	MACH3LOG_ERROR("Did not find histogram: {} in file: {}",ExpectedHistogramName,fDetectorSystematicsFileName);
	DetectorSystematicsFile->ls();
	throw MaCh3Exception(__FILE__, __LINE__);
      }
      RatioHistograms[iParam][iSample] = Histogram;
    }
  }

  //DB Now assign systematic weights
  for (size_t iEvent=0;iEvent<dunemcSamples.size();iEvent++) {
    int SampIndex = dunemcSamples[iEvent].SampleIndex;

    //DB Save each systematics ratio factor
    dunemcSamples[iEvent].DetectorSystematicRatios.resize(DetectorSystematicParameterNames.size());
    for (size_t iSyst=0;iSyst<DetectorSystematicParameterNames.size();iSyst++) {
      TH1* Histogram = RatioHistograms[iSyst][SampIndex];
      int HistogramBinIndex = -1;

      //DB Added support for 2D ratio histograms (and extension to 3D if we get that far). Current inputs are only in 1D so higher-D needs validating
      //   Assumes that 1D is reconstructed neutrino energy
      //   Assumes that 2D is reconstructed neutrino energy and direction
      if (Histogram->InheritsFrom(TH1::Class())) {
	HistogramBinIndex = Histogram->FindBin(dunemcSamples[iEvent].rw_erec);
      } else if (Histogram->InheritsFrom(TH2::Class())) {
	TH2* Histogram2D = static_cast<TH2*>(Histogram);
	HistogramBinIndex = Histogram2D->FindBin(dunemcSamples[iEvent].rw_erec,dunemcSamples[iEvent].rw_theta);
      } else {
	MACH3LOG_ERROR("Do not have support for 3D detector systematic binning yet");
	throw MaCh3Exception(__FILE__, __LINE__);
      }

      //DB Could probably do something smarter than just using ROOT histogram bin finding
      double Weight = RatioHistograms[iSyst][SampIndex]->GetBinContent(HistogramBinIndex);
      dunemcSamples[iEvent].DetectorSystematicRatios[iSyst] = Weight;
    }
  }
  
}

int SampleHandlerAtm::ReturnSampleIdentifier(std::vector<double> CVNScores, double MinDistanceToWall) {
  bool IsFullyContained = false;

  //DB Determine if the event is contained based on the config-defined option FCPCSeparation.
  //   If the minimum distance between any spacecharge hit and the TPC wall is larger than this value, it's fully contained
  if (MinDistanceToWall > 1e4 || MinDistanceToWall < 0) { //DB: ToDo Work out theoretical maximum
    return kEventSel_Unknown;
  } else if (MinDistanceToWall > FCPCSeparation) {
    IsFullyContained = true;
  } else {
    IsFullyContained = false;
  }

  //DB Current selections just defined whether it's e-like or mu-like
  //DB Not making the 0.55 and 0.56 magic numbers config-read because will eventually move to the argmax style selection (commented out on line above...)
  //int EventSelection = static_cast<int>(std::distance(CVNScores.begin(), max_element(CVNScores.begin(), CVNScores.end())));
  int EventSelection = kCVN_NC;
  if (CVNScores[kCVN_NuMu] > 0.56) { 
    EventSelection = kCVN_NuMu;
  } else if (CVNScores[kCVN_NuE] > 0.55) {
    EventSelection = kCVN_NuE;
  }

  //DB Use the above information and then return the index (which equates to the index of the sample in the SampleDetails)
  int SampIndex = kEventSel_Unknown;
  if (IsFullyContained) {
    if (EventSelection == kCVN_NuE)  {SampIndex = EventSelection_to_SampleIndex_Map[kEventSel_FC_NuE]; }
    if (EventSelection == kCVN_NuMu) {SampIndex = EventSelection_to_SampleIndex_Map[kEventSel_FC_NuMu];}
    if (EventSelection == kCVN_NC)   {SampIndex = EventSelection_to_SampleIndex_Map[kEventSel_FC_NC];  }    
  } else {
    if (EventSelection == kCVN_NuE)  {SampIndex = EventSelection_to_SampleIndex_Map[kEventSel_PC_NuE]; }
    if (EventSelection == kCVN_NuMu) {SampIndex = EventSelection_to_SampleIndex_Map[kEventSel_PC_NuMu];}
    if (EventSelection == kCVN_NC)   {SampIndex = EventSelection_to_SampleIndex_Map[kEventSel_PC_NC];  }    
  }
  
  return SampIndex;
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(const int KinPar, int iEvent) const {
  //DB For any new variable to plot, it needs to be added here
  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    return &(dunemcSamples[iEvent].enu_true);
  case kRecoNeutrinoEnergy:
    return &(dunemcSamples[iEvent].rw_erec);
  case kRecoHadronEnergy:
    return &(dunemcSamples[iEvent].rw_ehad);
  case kRecoLeptonEnergy:
    return &(dunemcSamples[iEvent].rw_elep);    
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
  //DB Use the above function just for ease
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}
