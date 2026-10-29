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

#include "TH2.h"
#include "TH3.h"

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
  //=============================================================
  //DB IsELike determines which reconstruction algorithm gets applied to which samples
  std::vector<std::string> EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);
  IsELike.resize(GetNSamples());
  BinningVariables.resize(GetNSamples());
  for(int iSample=0;iSample<GetNSamples();iSample++){
    const std::string TempTitle = EnabledSamples[iSample];
    IsELike[iSample] = Get<int>(SampleManager->raw()[TempTitle]["SampleOptions"]["IsELike"],__FILE__,__LINE__);

    std::vector<std::string> VarStr = (SampleManager->raw()[TempTitle]["Binning"]["VarStr"] && SampleManager->raw()[TempTitle]["Binning"]["VarStr"].IsScalar()) ?
      std::vector<std::string>{Get<std::string>(SampleManager->raw()[TempTitle]["Binning"]["VarStr"], __FILE__, __LINE__)} :
      Get<std::vector<std::string>>(SampleManager->raw()[TempTitle]["Binning"]["VarStr"],__FILE__, __LINE__);
    BinningVariables[iSample] = VarStr;
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
    std::string EigenInputFileMD5Sum = Get<std::string>(SampleManager->raw()["InputFiles"]["EigenFileMD5Sum"],__FILE__,__LINE__);

    if (EigenInputFile!="") {
      if (!CompareMD5Sum(EigenInputFileMD5Sum, EigenInputFile)) {
	MACH3LOG_ERROR("Eigen Input file MD5Sum does not match the expectation from the config file");
	throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  } else {
    EigenInputFile = "";
  }

  InputFileName = Get<std::string>(SampleManager->raw()["InputFiles"]["FileName"],__FILE__,__LINE__);
  InputFileDirectory = Get<std::string>(SampleManager->raw()["InputFiles"]["FileDirectory"],__FILE__,__LINE__);

  //=============================================================  
  // Per-event spline weights (optional)
  std::string InputSplineFile = Get<std::string>(SampleManager->raw()["InputFiles"]["InputSplines"],__FILE__,__LINE__);
  std::string SplineInputFileDirectory = Get<std::string>(SampleManager->raw()["InputFiles"]["InputSplineDirectory"],__FILE__,__LINE__);
  
  if (InputSplineFile!="") {
    std::string InputSplinesMD5Sum = Get<std::string>(SampleManager->raw()["InputFiles"]["InputSplinesMD5Sum"],__FILE__,__LINE__);

    if (InputSplineFile!="") {
      if (!CompareMD5Sum(InputSplinesMD5Sum, InputSplineFile)) {
	MACH3LOG_ERROR("Input Splines MD5Sum does not match the expectation from the config file");
	throw MaCh3Exception(__FILE__, __LINE__);
      }
    }

    fInputSplines.push_back(InputSplineFile);
  } else if (SplineInputFileDirectory!="") {
    //If reading multiple spline files
    std::filesystem::path directory = SplineInputFileDirectory;
    
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
      if (entry.is_regular_file() && entry.path().extension()==".root") {
	std::string FileName = entry.path();
	 fInputSplines.push_back(FileName);
      }
    }
  } else {
    MACH3LOG_INFO("No spline files provided");
  }

  //=============================================================  
  //DB Input file containing the detector systematic spectra ratios, with the expectation that HistogramTitle is formatted as SYSTEMATICNAME_SAMPLENAME
  if (SampleManager->raw()["InputFiles"]["DetectorSystematicsFileName"]) {
    fDetectorSystematicsFileName = Get<std::string>(SampleManager->raw()["InputFiles"]["DetectorSystematicsFileName"],__FILE__,__LINE__);
    std::string DetectorSystematicsFileMD5Sum = Get<std::string>(SampleManager->raw()["InputFiles"]["DetectorSystematicsFileMD5Sum"],__FILE__,__LINE__);
    if (fDetectorSystematicsFileName!="") {
      if (!CompareMD5Sum(DetectorSystematicsFileMD5Sum, fDetectorSystematicsFileName)) {
	MACH3LOG_ERROR("Detector Systematic Ratio file MD5Sum does not match the expectation from the config file");
	throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    
    //DB Read the expected dimension binning "configuration" from the config
    YAML::Node const& DetectorSystematicManager = (SampleManager->raw()["InputFiles"]["DetectorSystematicsConfiguration"]);
    if (!DetectorSystematicManager) {
      MACH3LOG_ERROR("[\"InputFiles\"][\"DetectorSystematicsConfiguration\"] not provided in Sample Config");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    for (const auto& item : DetectorSystematicManager) {
      std::string Key = item.first.as<std::string>();
      std::vector<int> KinematicParIndices;
      for (const auto& dimension : item.second) {
	std::string Dimension = dimension.as<std::string>();
	KinematicParIndices.push_back(KinematicParametersDUNE.at(Dimension));
      }
      DetectorSystematicDimensionMap[Key] = KinematicParIndices;
    }

    for (const auto& [key, values] : DetectorSystematicDimensionMap) {
      std::string tmpString = "";
      for (const auto& param : values) {
	tmpString += ReversedKinematicParametersDUNE.at(param)+", ";
      }
      MACH3LOG_INFO("Detector systematic {} configured with the following dimensions: {}",key,tmpString);
    }
  } else {
    fDetectorSystematicsFileName = "";
  }

  //=============================================================
  if (SampleManager->raw()["InputFiles"]["HONDAFluxUncerSplineFileName"]) {
    fHONDAFluxUncerSplineFileName = Get<std::string>(SampleManager->raw()["InputFiles"]["HONDAFluxUncerSplineFileName"],__FILE__,__LINE__);
    fHONDAFluxUncerSeparationValue = Get<double>(SampleManager->raw()["InputFiles"]["HONDAFluxUncerSeparationValue"],__FILE__,__LINE__);
    
    std::string HONDAFluxUncerSplineFileMD5Sum = Get<std::string>(SampleManager->raw()["InputFiles"]["HONDAFluxUncerSplineFileMD5Sum"],__FILE__,__LINE__);
    if (fHONDAFluxUncerSplineFileName!="") {
      if (!CompareMD5Sum(HONDAFluxUncerSplineFileMD5Sum, fHONDAFluxUncerSplineFileName)) {
	MACH3LOG_ERROR("HONDA Flux Uncertainty Spline file MD5Sum does not match the expectation from the config file");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }
  
  //=============================================================
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
    if (fInputSplines.size() == 0) {
      MACH3LOG_ERROR("Spline parameters defined in the config but no spline input files specified");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
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

  if (AtmFluxParameterNames.size()!=0) {
    RegisterIndividualFunctionalParameter(dunemcSamples, AtmFluxParameterNames, [](std::vector<double> const &, dunemc_atm &Event) {      
      Event.TotalAtmFluxSystematicWeight = 1.0;

      for (size_t iWeight=0;iWeight<Event.AtmFluxSystFactors.size();iWeight++) {
	Event.TotalAtmFluxSystematicWeight *= 1.0 + *Event.AtmFluxSystParVals[iWeight] * (Event.AtmFluxSystFactors[iWeight] - 1.0);
      }
    });
  }

}

void SampleHandlerAtm::AddAdditionalWeightPointers() {
  //DB ToDo: Try to "group" as many fixed weights as possible to reduce needless multiple multiplications
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].weights_flux));
    MCEvents[i].total_weight_pointers.push_back(&(ExposureScaling));
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].TotalDetectorSystematicWeight));
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].TotalAtmFluxSystematicWeight));
  }  
}

int SampleHandlerAtm::ReadFromEigen() {
  //DB Read the MC info from the Eigen binary blob -- but first check that the MD5 Checksum matches the expectation
  
  Eigen::MatrixXd Matrix;

  MACH3LOG_INFO("Reading Eigen matrix from:{}",EigenInputFile);
  ReadEigenMatrixFromFile(EigenInputFile.c_str(),Matrix);

  int nEntries = static_cast<int>(Matrix.rows());
  dunemcSamples.resize(nEntries);

  //DB The below needs modifying any time a new variable is read from, or derived from, the CAF file
  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    dunemcSamples[iEvent].neutrinoenergy_reco = Matrix(iEvent,Variables::neutrinoenergy_reco);
    dunemcSamples[iEvent].hadronicenergy_reco = Matrix(iEvent,Variables::hadronicenergy_reco);
    dunemcSamples[iEvent].leptonicenergy_reco = Matrix(iEvent,Variables::leptonicenergy_reco);
    dunemcSamples[iEvent].coszenith_reco = Matrix(iEvent,Variables::coszenith_reco);    
    dunemcSamples[iEvent].SampleIndex = static_cast<int>(Matrix(iEvent,Variables::SampleIndex));
    dunemcSamples[iEvent].neutrinoflavour_interacting = static_cast<int>(Matrix(iEvent,Variables::neutrinoflavour_interacting));
    dunemcSamples[iEvent].neutrinoflavour_generated = static_cast<int>(Matrix(iEvent,Variables::neutrinoflavour_generated));
    dunemcSamples[iEvent].OscChannelIndex = Matrix(iEvent,Variables::OscChannelIndex);
    dunemcSamples[iEvent].mode = Matrix(iEvent,Variables::mode);
    dunemcSamples[iEvent].IsCC = static_cast<int>(Matrix(iEvent,Variables::IsCC));
    dunemcSamples[iEvent].Target = static_cast<int>(Matrix(iEvent,Variables::Target));
    dunemcSamples[iEvent].neutrinoenergy_true = Matrix(iEvent,Variables::neutrinoenergy_true);
    dunemcSamples[iEvent].coszenith_true = Matrix(iEvent,Variables::coszenith_true);    
    dunemcSamples[iEvent].weights_flux = Matrix(iEvent,Variables::weights_flux);
    dunemcSamples[iEvent].MinDistToWall = Matrix(iEvent,Variables::MinDistToWall);
    dunemcSamples[iEvent].eventID = static_cast<uint>(Matrix(iEvent,Variables::eventID));
  }

  return nEntries;
}

void SampleHandlerAtm::TransferToEigen(std::string FileName) {
  Eigen::MatrixXd Matrix = Eigen::MatrixXd(dunemcSamples.size(),Variables::nVariables);

  //DB The below needs modifying any time a new variable is read from, or derived from, the CAF file
  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    Matrix(iEvent,Variables::neutrinoenergy_reco) = dunemcSamples[iEvent].neutrinoenergy_reco;
    Matrix(iEvent,Variables::hadronicenergy_reco) = dunemcSamples[iEvent].hadronicenergy_reco;
    Matrix(iEvent,Variables::leptonicenergy_reco) = dunemcSamples[iEvent].leptonicenergy_reco;
    Matrix(iEvent,Variables::coszenith_reco) = dunemcSamples[iEvent].coszenith_reco;
    Matrix(iEvent,Variables::SampleIndex) = dunemcSamples[iEvent].SampleIndex;
    Matrix(iEvent,Variables::neutrinoflavour_interacting) = dunemcSamples[iEvent].neutrinoflavour_interacting;
    Matrix(iEvent,Variables::neutrinoflavour_generated) = dunemcSamples[iEvent].neutrinoflavour_generated;
    Matrix(iEvent,Variables::OscChannelIndex) = dunemcSamples[iEvent].OscChannelIndex;
    Matrix(iEvent,Variables::mode) = dunemcSamples[iEvent].mode;
    Matrix(iEvent,Variables::IsCC) = dunemcSamples[iEvent].IsCC;    
    Matrix(iEvent,Variables::Target) = dunemcSamples[iEvent].Target;
    Matrix(iEvent,Variables::neutrinoenergy_true) = dunemcSamples[iEvent].neutrinoenergy_true;
    Matrix(iEvent,Variables::coszenith_true) = dunemcSamples[iEvent].coszenith_true;    
    Matrix(iEvent,Variables::weights_flux) = dunemcSamples[iEvent].weights_flux;
    Matrix(iEvent,Variables::MinDistToWall) = dunemcSamples[iEvent].MinDistToWall;
    Matrix(iEvent,Variables::eventID) = dunemcSamples[iEvent].eventID;        
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
      
      currentEvent_FromNuE.neutrinoenergy_reco = RecoENu;
      currentEvent_FromNuE.hadronicenergy_reco = RecoEHad;
      currentEvent_FromNuE.leptonicenergy_reco = RecoELep;            
      currentEvent_FromNuE.coszenith_reco = RecoCZ;
      currentEvent_FromNuE.SampleIndex = SampIndex;
      currentEvent_FromNuE.neutrinoflavour_interacting = InteractingPDG;
      currentEvent_FromNuE.neutrinoflavour_generated = (InteractingPDG > 0) ? 12 : -12;
      currentEvent_FromNuE.OscChannelIndex = static_cast<double>(GetOscChannel(OscillationChannels, currentEvent_FromNuE.neutrinoflavour_generated, currentEvent_FromNuE.neutrinoflavour_interacting));
      currentEvent_FromNuE.mode = M3Mode;
      currentEvent_FromNuE.IsCC = sr->mc.nu[0].iscc;
      currentEvent_FromNuE.Target = kTarget_Ar;
      currentEvent_FromNuE.neutrinoenergy_true = TrueNeutrinoEnergy;
      currentEvent_FromNuE.coszenith_true = -TrueNuMomentumVector.y(); // +Y in CAF files translates to +Z in typical CosZ
      currentEvent_FromNuE.weights_flux = xsec_w*flux_nue_w;
      currentEvent_FromNuE.MinDistToWall = MinDist;
      currentEvent_FromNuE.eventID = static_cast<uint>(iTreeEntry);
      
      struct dunemc_atm currentEvent_FromNuMu = currentEvent_FromNuE;
      
      currentEvent_FromNuMu.neutrinoflavour_generated = (InteractingPDG > 0) ? 14 : -14;
      currentEvent_FromNuMu.OscChannelIndex = static_cast<double>(GetOscChannel(OscillationChannels, currentEvent_FromNuMu.neutrinoflavour_generated, currentEvent_FromNuMu.neutrinoflavour_interacting));
      currentEvent_FromNuMu.weights_flux = xsec_w*flux_numu_w;

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
  SetupAtmFluxSplineSystematics();

  //DB Return the number of events loaded
  return static_cast<int>(dunemcSamples.size());
}

void SampleHandlerAtm::SetupMC() {
  //DB Pass the Core-required information over 
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent) {
    MCEvents[iEvent].enu_true = dunemcSamples[iEvent].neutrinoenergy_true;
    MCEvents[iEvent].isNC = !dunemcSamples[iEvent].IsCC;
    MCEvents[iEvent].nupdg = dunemcSamples[iEvent].neutrinoflavour_interacting;
    MCEvents[iEvent].nupdgUnosc = dunemcSamples[iEvent].neutrinoflavour_generated;
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

      //DB Check that the dimensions for this parameter are in the configured map
      if (DetectorSystematicDimensionMap.count(ParamName) == 0) {
	MACH3LOG_ERROR("Detector Systematic Parameter {} has no configuration to determine what the dimensions of the ratio spectra are binned in",ParamName);
	throw MaCh3Exception(__FILE__, __LINE__);
      }
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

      //DB Check how many dimensions were configured
      std::vector<int> ConfiguredDimensions = DetectorSystematicDimensionMap.at(DetectorSystematicParameterNames[iSyst]);
      size_t nConfiguredDimensions = ConfiguredDimensions.size();

      //DB Added support for higher-dimension ratio histograms (and extension to 3D if we get that far). Current inputs are only in 1D so higher-D needs validating
      //   Assumes that 1D is reconstructed neutrino energy
      //   Assumes that 2D is reconstructed neutrino energy and direction
      //   Assumes that 3D is reconstructed hadronic energy, leptonic energy, and direction      
      if (Histogram->GetDimension() == 1) {
	if (nConfiguredDimensions != 1) {
	  MACH3LOG_ERROR("Misconfigured Detector Systematics Parameter {} - Config has given {} dimension for the spectra histogram, whereas the histogram itself has 1 dimension",DetectorSystematicParameterNames[iSyst],nConfiguredDimensions);
	  throw MaCh3Exception(__FILE__, __LINE__);
	}
	HistogramBinIndex = Histogram->FindBin(*GetPointerToKinematicParameter(ConfiguredDimensions[0],static_cast<int>(iEvent)));
      } else if (Histogram->GetDimension() == 2) {
	if (nConfiguredDimensions != 2) {
          MACH3LOG_ERROR("Misconfigured Detector Systematics Parameter {} - Config has given {} dimension for the spectra histogram, whereas the histogram itself has 2 dimensions",DetectorSystematicParameterNames[iSyst],nConfiguredDimensions);
          throw MaCh3Exception(__FILE__, __LINE__);
	}
	TH2* Histogram2D = static_cast<TH2*>(Histogram);
	HistogramBinIndex = Histogram2D->FindBin(*GetPointerToKinematicParameter(ConfiguredDimensions[0],static_cast<int>(iEvent)),*GetPointerToKinematicParameter(ConfiguredDimensions[1],static_cast<int>(iEvent)));
      } else if (Histogram->GetDimension() == 3) {
	if (nConfiguredDimensions != 3) {
          MACH3LOG_ERROR("Misconfigured Detector Systematics Parameter {} - Config has given {} dimension for the spectra histogram, whereas the histogram itself has 3 dimensions",DetectorSystematicParameterNames[iSyst],nConfiguredDimensions);
          throw MaCh3Exception(__FILE__, __LINE__);
	}
	TH3* Histogram3D = static_cast<TH3*>(Histogram);
	HistogramBinIndex = Histogram3D->FindBin(*GetPointerToKinematicParameter(ConfiguredDimensions[0],static_cast<int>(iEvent)),*GetPointerToKinematicParameter(ConfiguredDimensions[1],static_cast<int>(iEvent)),*GetPointerToKinematicParameter(ConfiguredDimensions[2],static_cast<int>(iEvent)));
      } else {
	MACH3LOG_ERROR("Found detector systematic ratio histogram which has more than 3 dimensions");
	throw MaCh3Exception(__FILE__, __LINE__);
      }

      //DB Could probably do something smarter than just using ROOT histogram bin finding
      double Weight = RatioHistograms[iSyst][SampIndex]->GetBinContent(HistogramBinIndex);
      dunemcSamples[iEvent].DetectorSystematicRatios[iSyst] = Weight;
    }
  }
  
}

void SampleHandlerAtm::SetupAtmFluxSplineSystematics() {
  //DB Firstly set detector systematic weight to 1.0
  for (size_t iEvent=0;iEvent<dunemcSamples.size();iEvent++) {
    dunemcSamples[iEvent].TotalAtmFluxSystematicWeight = 1.0;
  }

  //DB Figure out whether we have detector systematics configured -- required that any Detector Systematic parameter's group is "AtmFlux"                                                                  
  int NParams = ParHandler->GetNParameters();
  for (int iParam=0;iParam<NParams;iParam++) {
    if (ParHandler->IsParFromGroup(iParam,"AtmFlux")) {
      std::string ParamName = ParHandler->GetParFancyName(iParam);
      AtmFluxParameterNames.push_back(ParamName);
    }
  }

  if (AtmFluxParameterNames.size()==0) {
    return;
  }

  //====================================================================
  
  bool HasHondaFluxUncertainty =
    std::find(AtmFluxParameterNames.begin(), AtmFluxParameterNames.end(),
              "HondaFlux_TotalUncertainty_Below1GeV") != AtmFluxParameterNames.end() &&
    std::find(AtmFluxParameterNames.begin(), AtmFluxParameterNames.end(),
              "HondaFlux_TotalUncertainty_Above1GeV") != AtmFluxParameterNames.end();
  
  if ((HasHondaFluxUncertainty && fHONDAFluxUncerSplineFileName=="") || (!HasHondaFluxUncertainty && fHONDAFluxUncerSplineFileName!="")) {
    MACH3LOG_ERROR("Bad config");
    throw MaCh3Exception(__FILE__, __LINE__);      
  }

  if (fHONDAFluxUncerSplineFileName!="") {
    TFile* File = new TFile(fHONDAFluxUncerSplineFileName.c_str());
    if (!File || File->IsZombie()) {
      MACH3LOG_ERROR("Did not find Honda Flux Uncertainty Spline file: {}",fHONDAFluxUncerSplineFileName);
      throw MaCh3Exception(__FILE__, __LINE__);      
    }

    TSpline3* LowEnergySpline = File->Get<TSpline3>("HondaTotalUncertainty_EnergyBin_0");
    if (!LowEnergySpline) {
      MACH3LOG_ERROR("Did not find Honda Flux Uncertainty Spline: {}","HondaTotalUncertainty_EnergyBin_0");
      throw MaCh3Exception(__FILE__, __LINE__);            
    }
    TSpline3* HighEnergySpline = File->Get<TSpline3>("HondaTotalUncertainty_EnergyBin_1");
    if (!HighEnergySpline) {
      MACH3LOG_ERROR("Did not find Honda Flux Uncertainty Spline: {}","HondaTotalUncertainty_EnergyBin_1");
      throw MaCh3Exception(__FILE__, __LINE__);            
    }

    int LowEnergySystParIndex = -1;
    int HighEnergySystParIndex = -1;

    for (int iParam=0;iParam<NParams;iParam++) {
      if (ParHandler->GetParFancyName(iParam) == "HondaFlux_TotalUncertainty_Below1GeV") LowEnergySystParIndex = iParam;
      if (ParHandler->GetParFancyName(iParam) == "HondaFlux_TotalUncertainty_Above1GeV") HighEnergySystParIndex = iParam;      
    }

    for (size_t iEvent=0;iEvent<dunemcSamples.size();iEvent++) {
      dunemcSamples[iEvent].AtmFluxSystFactors.resize(1);
      dunemcSamples[iEvent].AtmFluxSystParVals.resize(1);      
      
      double Factor = -1;
      int ParIndex = -1;
      
      if (dunemcSamples[iEvent].neutrinoenergy_true < fHONDAFluxUncerSeparationValue) {
	Factor = LowEnergySpline->Eval(dunemcSamples[iEvent].neutrinoenergy_true);
	ParIndex = LowEnergySystParIndex;
      } else {
	Factor = HighEnergySpline->Eval(dunemcSamples[iEvent].neutrinoenergy_true);
	ParIndex = HighEnergySystParIndex;	
      }
      
      dunemcSamples[iEvent].AtmFluxSystFactors[0] = Factor;
      dunemcSamples[iEvent].AtmFluxSystParVals[0] = ParHandler->RetPointer(ParIndex);
    }
  }

  //====================================================================    
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
    return &(dunemcSamples[iEvent].neutrinoenergy_true);
  case kRecoNeutrinoEnergy:
    return &(dunemcSamples[iEvent].neutrinoenergy_reco);
  case kRecoHadronEnergy:
    return &(dunemcSamples[iEvent].hadronicenergy_reco);
  case kRecoLeptonEnergy:
    return &(dunemcSamples[iEvent].leptonicenergy_reco);    
  case kTrueCosZ:
    return &(dunemcSamples[iEvent].coszenith_true);
  case kRecoCosZ:
    return &(dunemcSamples[iEvent].coszenith_reco);
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
