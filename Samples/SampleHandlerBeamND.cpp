#include "SampleHandlerBeamND.h"

//Here nullptr is passed instead of OscCov to prevent oscillation calculations from being performed for the ND Samples
SampleHandlerBeamND::SampleHandlerBeamND(std::string mc_version_, ParameterHandlerGeneric* ParHandler_,  TMatrixD* nd_cov_) : SampleHandlerFD(mc_version_, ParHandler_) {
  if(!nd_cov_){
    MACH3LOG_ERROR("You've passed me a nullptr to a ND covarince matrix... ");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  NDCovMatrix = nd_cov_;

  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerBeamND::~SampleHandlerBeamND() {
}

void SampleHandlerBeamND::Init() {
  dunendmcSamples.resize(nSamples,dunemc_beamnd());
  
  IsFHC = SampleManager->raw()["DUNESampleBools"]["isFHC"].as<double>();
  iselike = SampleManager->raw()["DUNESampleBools"]["iselike"].as<bool>();
  pot = SampleManager->raw()["POT"].as<double>();

  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void SampleHandlerBeamND::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(new BinnedSplineHandlerDUNE(ParHandler,Modes.get()));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

void SampleHandlerBeamND::SetupWeightPointers() {
  for (size_t i = 0; i < dunendmcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.push_back(&(dunendmcSamples[i].pot_s));
    MCSamples[i].total_weight_pointers.push_back(&(dunendmcSamples[i].norm_s));
    MCSamples[i].total_weight_pointers.push_back(&(dunendmcSamples[i].rw_berpaacvwgt));
    MCSamples[i].total_weight_pointers.push_back(&(dunendmcSamples[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));
  }
}

int SampleHandlerBeamND::SetupExperimentMC() {

  // dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  
  MACH3LOG_INFO("-------------------------------------------------------------------");

  if (IsFHC) { 
    norm_s = (1e21/1.905e21);
  } else {
    norm_s = (1e21/1.5e21);
  }
  pot_s = (pot)/1e21;

  TChain* _data = new TChain("caf");
  for (size_t iSample=0;iSample<mc_files.size();iSample++) {
    MACH3LOG_INFO("Adding file to TChain: {}", mc_files[iSample]);
    // HH: currently storing the norm_s and POT_s in a map since this is a sample/file specific thing, maybe this could be done in a more elegant way
    norm_map[mc_files[iSample]] = std::vector<double>{norm_s, pot_s};

    // HH: Check whether the file exists, see https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
    int _add_rtn = _data->Add(mc_files[iSample].c_str(), -1);
    if(_add_rtn == 0){
      MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", mc_files[iSample]);
      throw MaCh3Exception(__FILE__, __LINE__);
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
  dunendmcSamples.resize(nEntries);
  _data->GetEntry(0);

  //FILL DUNE STRUCT
  for (unsigned int i = 0; i < nEntries; ++i) { // Loop through tree
    if (i % (nEntries/10) == 0) {
      MACH3LOG_INFO("Processing entry {}/{} ({}%)", i, nEntries, (i*100)/nEntries);
    }
    _data->GetEntry(i);

    std::string CurrFileName = _data->GetCurrentFile()->GetName();
    dunendmcSamples[i].nupdgUnosc = _nuPDGunosc;
    dunendmcSamples[i].nupdg = _nuPDG;
    dunendmcSamples[i].OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, dunendmcSamples[i].nupdgUnosc, dunendmcSamples[i].nupdg));

    // POT stuff
    dunendmcSamples[i].norm_s = norm_s;
    dunendmcSamples[i].pot_s = pot_s;
    
    dunendmcSamples[i].rw_erec = _erec;
    dunendmcSamples[i].rw_erec_shifted = _erec;
    dunendmcSamples[i].rw_erec_lep = _erec_lep;
    dunendmcSamples[i].rw_erec_had = (_erec - _erec_lep);
    dunendmcSamples[i].rw_yrec = ((_erec - _erec_lep)/_erec);
    dunendmcSamples[i].rw_etru = _ev; // in GeV
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


const double* SampleHandlerBeamND::GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent) {
  double* KinematicValue;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunendmcSamples[iEvent].rw_etru);
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

const double* SampleHandlerBeamND::GetPointerToKinematicParameter(double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

const double* SampleHandlerBeamND::GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

double SampleHandlerBeamND::ReturnKinematicParameter(int KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamND::ReturnKinematicParameter(std::string KinematicParameter, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iEvent);
}

void SampleHandlerBeamND::SetupFDMC() {
  // dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  // FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);
  
  for(int iEvent = 0 ;iEvent < int(GetNEvents()); ++iEvent){
    MCSamples[iEvent].rw_etru = &(dunendmcSamples[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunendmcSamples[iEvent].mode);
    MCSamples[iEvent].Target = &(dunendmcSamples[iEvent].Target); 
    MCSamples[iEvent].isNC = !(dunendmcSamples[iEvent].rw_isCC);
    MCSamples[iEvent].nupdgUnosc = &(dunendmcSamples[iEvent].nupdgUnosc);
    MCSamples[iEvent].nupdg = &(dunendmcSamples[iEvent].nupdg);
  }
}

// Set the covariance matrix for this class
void SampleHandlerBeamND::setNDCovMatrix() {
  if (NDCovMatrix == NULL) {
    std::cerr << "Could not find ND Detector covariance matrix you provided to setCovMatrix" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  int nXBinsNDCov = static_cast<int>(Binning.XBinEdges.size()-1);
  int nYBinsNDCov = static_cast<int>(Binning.YBinEdges.size()-1);
  int covSize = nXBinsNDCov*nYBinsNDCov;

  if (covSize != NDCovMatrix->GetNrows()) {
    std::cerr << "Sample dimensions do not match ND Detector Covariance!" << std::endl;
    std::cerr << "Sample XBinsNDCov * YBinsNDCov = " << covSize  << std::endl;
    std::cerr << "ND Detector Covariance = " << NDCovMatrix->GetNrows() << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  
  std::vector<double> FlatCV;
  int iter = 0;

  // 2D -> 1D Array
  for (int xBin = 0; xBin < nXBinsNDCov; xBin++) 
  {
    for (int yBin = 0; yBin < nYBinsNDCov; yBin++) 
    {
        const int idx = Binning.GetBinSafe(xBin, yBin);
        double CV = SampleHandlerFD_data[idx];
        FlatCV.push_back(CV);

        if(CV>0) (*NDCovMatrix)(iter,iter) += 1/CV;

	    iter++;
	}
  }

  NDInvertCovMatrix = new double*[covSize]();
  // Set the defaults to true
  for(int i = 0; i < covSize; i++)
  {
    NDInvertCovMatrix[i] = new double[covSize]();
    for (int j = 0; j < covSize; j++)
    {
        NDInvertCovMatrix[i][j] = 0.;
    }
  }

  TMatrixD* NDInvCovMatrix=static_cast<TMatrixD*>(NDCovMatrix->Clone());
  NDInvCovMatrix->Invert();

 
  //Scale back to inverse absolute cov and use standard double
  for (int i = 0; i < covSize; i++) {
    for (int j = 0; j < covSize; ++j) {
      const double f = FlatCV[i] * FlatCV[j];
      if(f != 0) NDInvertCovMatrix[i][j] = (*NDInvCovMatrix)(i,j)/f;
      else NDInvertCovMatrix[i][j] = 0.;
	}
  }

}

//New likelihood calculation for ND samples using detector covariance matrix
double SampleHandlerBeamND::GetLikelihood()
{
  if (SampleHandlerFD_data == NULL) {
      std::cerr << "data sample is empty!" << std::endl;
      return -1;
  }

  //if (NDInvertCovMatrix == NULL) {
  //setNDCovMatrix();
  //}

  if (!isNDCovSet) {
    setNDCovMatrix();
    isNDCovSet = true;
  }

  int nXBinsNDCov = static_cast<int>(Binning.XBinEdges.size()-1);
  int nYBinsNDCov = static_cast<int>(Binning.YBinEdges.size()-1);

  int covSize = nXBinsNDCov*nYBinsNDCov;

  std::vector<double> FlatData;
  std::vector<double> FlatMCPred;

  //2D -> 1D 
  for (int xBin = 0; xBin < nXBinsNDCov; xBin++) 
  {
    for (int yBin = 0; yBin < nYBinsNDCov; yBin++) 
    {
      const int idx = Binning.GetBinSafe(xBin, yBin);
      double MCPred = SampleHandlerFD_array[idx];
      FlatMCPred.push_back(MCPred);

      double DataVal = SampleHandlerFD_data[idx];
      FlatData.push_back(DataVal);
    }
  }


  double negLogL = 0.;
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:negLogL)
#endif

  for (int i = 0; i < covSize; i++) {
    for (int j = 0; j <= i; ++j) {
        //KS: Since matrix is symetric we can calcaute non daigonal elements only once and multiply by 2, can bring up to factor speed decrease.   
        int scale = 1;
        if(i != j) scale = 2;
        negLogL += scale * 0.5*(FlatData[i] - FlatMCPred[i])*(FlatData[j] - FlatMCPred[j])*NDInvertCovMatrix[i][j];
      }
  }

  return negLogL;
}
