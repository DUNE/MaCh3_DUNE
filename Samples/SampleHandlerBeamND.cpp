#include "SampleHandlerBeamND.h"
#include "Samples/StructsDUNE.h"

//Here nullptr is passed instead of OscCov to prevent oscillation calculations from being performed for the ND Samples
SampleHandlerBeamND::SampleHandlerBeamND(std::string mc_version_, ParameterHandlerGeneric* ParHandler_, BeamNDCov beamNDCov_) : SampleHandlerBase(mc_version_, ParHandler_) {
  if(!(beamNDCov_.NDCov_FHC && beamNDCov_.NDCov_RHC && beamNDCov_.NDCov_all)){
    MACH3LOG_ERROR("You've passed me a nullptr to a ND covarince matrix... ");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  beamNDCov = beamNDCov_;

  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerBeamND::~SampleHandlerBeamND() {
}

void SampleHandlerBeamND::Init() {
  beamNDSampleDetails.resize(GetNSamples());
  
  auto EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);

  for (int i = 0; i < GetNSamples(); i++){
    const auto TempTitle = EnabledSamples[i];
    beamNDSampleDetails[i].isFHC = SampleManager->raw()[TempTitle]["DUNESampleBools"]["isFHC"].as<double>();
    beamNDSampleDetails[i].iselike = SampleManager->raw()[TempTitle]["DUNESampleBools"]["iselike"].as<bool>();
    beamNDSampleDetails[i].pot = SampleManager->raw()[TempTitle]["POT"].as<double>();

    if (beamNDSampleDetails[i].isFHC) { 
      beamNDSampleDetails[i].norm_s = (1e21/1.905e21);
    } else {
      beamNDSampleDetails[i].norm_s = (1e21/1.5e21);
    }
    beamNDSampleDetails[i].pot_s = (beamNDSampleDetails[i].pot)/1e21;

    MACH3LOG_INFO("Setting up beam ND sample {}", GetSampleTitle(i));
    MACH3LOG_INFO("- isFHC: {}", beamNDSampleDetails[i].isFHC);
    MACH3LOG_INFO("- iselike: {}", beamNDSampleDetails[i].iselike);
  }

  MACH3LOG_INFO("-------------------------------------------------------------------");
}

// ************************************************
void SampleHandlerBeamND::InititialiseData()
{
  // ************************************************
  // Reweight MC to match
  Reweight();
  // set asimov data
  for (int iSample = 0; iSample < GetNSamples(); iSample++)
  {
    AddData(iSample, GetMCArray(iSample));
  }
}

void SampleHandlerBeamND::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(new BinnedSplineHandlerDUNE(ParHandler,Modes.get()));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

void SampleHandlerBeamND::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < dunendmcSamples.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(dunendmcSamples[i].pot_s));
    MCEvents[i].total_weight_pointers.push_back(&(dunendmcSamples[i].norm_s));
    MCEvents[i].total_weight_pointers.push_back(&(dunendmcSamples[i].rw_berpaacvwgt));
    MCEvents[i].total_weight_pointers.push_back(&(dunendmcSamples[i].flux_w));
  }
}

int SampleHandlerBeamND::SetupExperimentMC() {

  // dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  
  MACH3LOG_INFO("-------------------------------------------------------------------");

  TChain* _data = new TChain("caf");
  // Maps the file index within the TChain (GetTreeNumber()) to its sample index.
  std::vector<size_t> fileIndexToSample;
  for (size_t iSample=0;iSample<SampleDetails.size();iSample++) {
    for (const std::vector<std::string>& files : SampleDetails[iSample].mc_files) {
      for (const std::string& filename : files) {

        MACH3LOG_INFO("Adding file to TChain: {}", filename);
        // HH: Check whether the file exists, see https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
        int _add_rtn = _data->Add(filename.c_str(), -1);
        if(_add_rtn == 0){
          MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", filename);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        // Each file added (glob patterns may expand to multiple) gets the same sample index.
        // We query how many files are now in the chain to know how many entries to push.
        const int nFilesNow = _data->GetListOfFiles()->GetEntries();
        while(static_cast<int>(fileIndexToSample.size()) < nFilesNow){
          fileIndexToSample.push_back(iSample);
        }
      }
    }
  }

  if(_data){
    MACH3LOG_INFO("Number of entries in TChain: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Failed to create the TChain.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
  _data->SetBranchStatus("Ev_reco", 1);
  _data->SetBranchAddress("Ev_reco", &_erec);
  _data->SetBranchStatus("Elep_reco", 1);
  _data->SetBranchAddress("Elep_reco", &_erec_lep);
  _data->SetBranchStatus("mode",1);
  _data->SetBranchAddress("mode",&_mode);
  _data->SetBranchStatus("isCC", 1);
  _data->SetBranchAddress("isCC", &_isCC);
  _data->SetBranchStatus("reco_q", 1);
  _data->SetBranchAddress("reco_q", &_reco_q);
  _data->SetBranchStatus("nuPDGunosc", 1);
  _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
  _data->SetBranchStatus("nuPDG", 1);
  _data->SetBranchAddress("nuPDG", &_nuPDG);
  _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
  _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);

  size_t nEntries = static_cast<size_t>(_data->GetEntries());
  size_t countwidth = nEntries / 10;
  dunendmcSamples.resize(nEntries);
  _data->GetEntry(0);

  //FILL DUNE STRUCT
  for (unsigned int i = 0; i < nEntries; ++i) { // Loop through tree
    _data->GetEntry(i);

    if (i % countwidth == 0) {
      M3::Utils::PrintProgressBar(i, static_cast<Long64_t>(nEntries));
    }

    const Int_t treeNum = _data->GetTreeNumber();
    if(treeNum < 0 || static_cast<size_t>(treeNum) >= fileIndexToSample.size()){
      MACH3LOG_ERROR("GetTreeNumber() returned {} which is out of range [0, {})", treeNum, fileIndexToSample.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const size_t sample_index = fileIndexToSample[static_cast<size_t>(treeNum)];
    dunendmcSamples[i].SampleIndex = static_cast<unsigned int>(sample_index);
    dunendmcSamples[i].nupdgUnosc = _nuPDGunosc;
    dunendmcSamples[i].nupdg = _nuPDG;
    dunendmcSamples[i].OscChannelIndex = static_cast<double>(GetOscChannel(SampleDetails[sample_index].OscChannels, dunendmcSamples[i].nupdgUnosc, dunendmcSamples[i].nupdg));

    // POT stuff
    dunendmcSamples[i].norm_s = beamNDSampleDetails[sample_index].norm_s;
    dunendmcSamples[i].pot_s = beamNDSampleDetails[sample_index].pot_s;
    
    dunendmcSamples[i].rw_erec = _erec;
    dunendmcSamples[i].rw_erec_shifted = _erec;
    dunendmcSamples[i].rw_erec_lep = _erec_lep;
    dunendmcSamples[i].rw_erec_had = (_erec - _erec_lep);
    dunendmcSamples[i].rw_yrec = ((_erec - _erec_lep)/_erec);
    dunendmcSamples[i].enu_true = _ev; // in GeV
    dunendmcSamples[i].rw_theta = _LepNuAngle;
    dunendmcSamples[i].rw_isCC = _isCC;
    dunendmcSamples[i].rw_reco_q = _reco_q;
    dunendmcSamples[i].rw_berpaacvwgt = _BeRPA_cvwgt;
    
    //Assume everything is on Argon for now....
    dunendmcSamples[i].Target = 40;

    int M3Mode = Modes->GetModeFromGenerator(std::abs(_mode));
    if (!_isCC) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    dunendmcSamples[i].mode = M3Mode;
    
    dunendmcSamples[i].flux_w = 1.0;
  }
  
  //_sampleFile->Close();
  _data->Reset();
  delete _data;
  return static_cast<int>(nEntries);
  
}


const double* SampleHandlerBeamND::GetPointerToKinematicParameter(const int KinPar, const int iEvent) const{
  const double* KinematicValue;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunendmcSamples[iEvent].enu_true);
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &(dunendmcSamples[iEvent].rw_erec_shifted);
    break;
  case kyRec:
    KinematicValue = &(dunendmcSamples[iEvent].rw_yrec);
    break;
  case kOscChannel:
    KinematicValue = &(dunendmcSamples[iEvent].OscChannelIndex);
    break;
  case kMode:
    KinematicValue = &(dunendmcSamples[iEvent].mode);
    break;
  case kIsFHC:
    KinematicValue = &(IsFHC);
    break;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}



double SampleHandlerBeamND::ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

void SampleHandlerBeamND::SetupMC() {
  // dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  // FarDetectorCoreInfo *fdobj = &(MCEvents[iSample]);
  
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    MCEvents[iEvent].enu_true = dunendmcSamples[iEvent].enu_true;
    MCEvents[iEvent].isNC = !(dunendmcSamples[iEvent].rw_isCC);
    MCEvents[iEvent].nupdgUnosc = dunendmcSamples[iEvent].nupdgUnosc;
    MCEvents[iEvent].nupdg = dunendmcSamples[iEvent].nupdg;
    MCEvents[iEvent].NominalSample = dunendmcSamples[iEvent].SampleIndex;
  }
}

// Set the covariance matrix for this class
void SampleHandlerBeamND::setNDCovMatrix() const {
  const int covSize = Binning->GetNBins();

  // Build a working covSize x covSize matrix:
  // if useCombinedNDCov: copy NDCov_all directly (spans all sub-samples)
  // else: assemble a block-diagonal from NDCov_FHC + NDCov_RHC
  TMatrixD WorkCov(covSize, covSize);
  WorkCov.Zero();

  if (beamNDCov.useCombinedNDCov) {
    if (beamNDCov.NDCov_all == nullptr) {
      MACH3LOG_ERROR("NDCov_all is nullptr");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (covSize != beamNDCov.NDCov_all->GetNrows()) {
      MACH3LOG_ERROR("Combined ND covariance size {} does not match total bins {}",
                     beamNDCov.NDCov_all->GetNrows(), covSize);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    WorkCov = *(beamNDCov.NDCov_all);
  }

  std::vector<double> FlatCV(covSize);
  int globalBin = 0;

  /*==========================================================================================================================
  //DB Do not trust the below code (between the lines of '=') - changed just to compile  
  */
  
  for (int iSample = 0; iSample < static_cast<int>(SampleDetails.size()); iSample++) {
    const int nBins = Binning->GetNBins(iSample);
    const int blockSize = nBins;

    if (!beamNDCov.useCombinedNDCov) {
      // Fill the block-diagonal entry for this sub-sample
      const bool isFHCSample = (beamNDSampleDetails[iSample].isFHC != 0);
      TMatrixD* SampleCov = isFHCSample ? beamNDCov.NDCov_FHC : beamNDCov.NDCov_RHC;
      if (SampleCov == nullptr) {
        MACH3LOG_ERROR("NDCov matrix for sub-sample {} is nullptr", iSample);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      if (blockSize != SampleCov->GetNrows()) {
        MACH3LOG_ERROR("Covariance size {} does not match block size {} for sub-sample {}",
                       SampleCov->GetNrows(), blockSize, iSample);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      for (int i = 0; i < blockSize; i++)
        for (int j = 0; j < blockSize; j++)
          WorkCov(globalBin + i, globalBin + j) = (*SampleCov)(i, j);
    }

    // Fill flat CV vector and add statistical term to diagonal
    int localBin = 0;
    for (int iBin = 0; iBin < nBins; iBin++) {
      const double CV = SampleHandler_data[iBin];
      FlatCV[globalBin + localBin] = CV;
      if (CV > 0)
	WorkCov(globalBin + localBin, globalBin + localBin) += 1.0 / CV;
      localBin++;
    }
    globalBin += blockSize;
  }

  //==========================================================================================================================

  // Invert the working matrix
  WorkCov.Invert();

  // Allocate and fill NDInvertCovMatrix, scaling back to absolute inverse cov
  NDInvertCovMatrix = new double *[covSize]();
  for (int i = 0; i < covSize; i++) {
    NDInvertCovMatrix[i] = new double[covSize]();
    for (int j = 0; j < covSize; j++) {
      const double f = FlatCV[i] * FlatCV[j];
      NDInvertCovMatrix[i][j] = (f != 0) ? WorkCov(i, j) / f : 0.;
    }
  }
}

// New likelihood calculation for ND samples using detector covariance matrix
double SampleHandlerBeamND::GetLikelihood() const {
  if (SampleHandler_data.empty()) {
    MACH3LOG_ERROR("data sample is empty!");
    return -1;
  }

  if (!isNDCovSet) {
    setNDCovMatrix();
    isNDCovSet = true;
  }

  const int covSize = Binning->GetNBins();

  std::vector<double> FlatData(covSize);
  std::vector<double> FlatMCPred(covSize);

  /*==========================================================================================================================
  //DB Do not trust the below code (between the lines of '=') - changed just to compile
  */
  
  // 2D -> 1D, iterating over all sub-samples in the same order as setNDCovMatrix
  for (int iSample = 0; iSample < static_cast<int>(SampleDetails.size()); iSample++) {
    for (int iBin = 0; iBin < Binning->GetNBins(iSample); iBin++) {
      FlatData[iBin] = SampleHandler_data[iBin];
      FlatMCPred[iBin] = SampleHandler_array[iBin];
    }
  }

  //==========================================================================================================================

  double negLogL = 0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+ : negLogL)
#endif

  for (int i = 0; i < covSize; i++) {
    for (int j = 0; j <= i; ++j) {
      const int scale = (i != j) ? 2 : 1;
      negLogL += scale * 0.5 * (FlatData[i] - FlatMCPred[i]) *
                 (FlatData[j] - FlatMCPred[j]) * NDInvertCovMatrix[i][j];
    }
  }

  return negLogL;
}
