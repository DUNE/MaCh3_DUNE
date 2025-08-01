#include "SampleHandlerBeamND.h"

//Here nullptr is passed instead of OscCov to prevent oscillation calculations from being performed for the ND Samples
SampleHandlerBeamND::SampleHandlerBeamND(std::string mc_version_, ParameterHandlerGeneric* ParHandler_,  TMatrixD* nd_cov_, ParameterHandlerOsc* osc_cov_=nullptr) : SampleHandlerFD(mc_version_, ParHandler_, osc_cov_) {
  OscParHandler = nullptr;
  
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
  dunendmcSamples.resize(nSamples,dunemc_base());
  
  IsFHC = SampleManager->raw()["DUNESampleBools"]["isFHC"].as<double>();
  iselike = SampleManager->raw()["DUNESampleBools"]["iselike"].as<bool>();
  if (CheckNodeExists(SampleManager->raw(), "Exposure")) {
    double exposure = SampleManager->raw()["Exposure"].as<double>();
    if (CheckNodeExists(SampleManager->raw(), "ExposureToPOT")) {
      double exposure_to_pot = SampleManager->raw()["ExposureToPOT"].as<double>();
      pot = exposure * exposure_to_pot;
    } else {
      MACH3LOG_ERROR("ExposureToPOT not defined in {}, please add this!", SampleManager->GetFileName());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    MACH3LOG_INFO("Using {} MWyr exposure, which gives {} POT", exposure, pot);
  } else if (CheckNodeExists(SampleManager->raw(), "POT")) {
    pot = SampleManager->raw()["POT"].as<double>();
    MACH3LOG_INFO("Using {} POT rather than exposure", pot);
  } else {
    MACH3LOG_ERROR("Exposure or POT not defined in {}, please add this!", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void SampleHandlerBeamND::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(new BinnedSplineHandlerDUNE(ParHandler,Modes));
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
    for (int j = 0; j < dunendmcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 5;
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunendmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = &(dunendmcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][3] = &(dunendmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

int SampleHandlerBeamND::SetupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files.at(iSample));
  
  TChain* _data = new TChain("caf");
  _data->Add(mc_files.at(iSample).c_str());

  if(_data){
    MACH3LOG_INFO("Found \"caf\" tree in {}", mc_files[iSample]);
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find \"caf\" tree in {}", mc_files[iSample]);
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

  _data->SetBranchStatus("eRecoP", 1);
  _data->SetBranchAddress("eRecoP", &_eRecoP);
  _data->SetBranchStatus("eRecoPip", 1);
  _data->SetBranchAddress("eRecoPip", &_eRecoPip);
  _data->SetBranchStatus("eRecoPim", 1);
  _data->SetBranchAddress("eRecoPim", &_eRecoPim);
  _data->SetBranchStatus("eRecoPi0", 1);
  _data->SetBranchAddress("eRecoPi0", &_eRecoPi0);
  _data->SetBranchStatus("eRecoN", 1);
  _data->SetBranchAddress("eRecoN", &_eRecoN);

  _data->SetBranchStatus("LepE", 1);
  _data->SetBranchAddress("LepE", &_LepE);
  _data->SetBranchStatus("eP", 1);
  _data->SetBranchAddress("eP", &_eP);
  _data->SetBranchStatus("ePip", 1);
  _data->SetBranchAddress("ePip", &_ePip);
  _data->SetBranchStatus("ePim", 1);
  _data->SetBranchAddress("ePim", &_ePim);
  _data->SetBranchStatus("ePi0", 1);
  _data->SetBranchAddress("ePi0", &_ePi0);
  _data->SetBranchStatus("eN", 1);
  _data->SetBranchAddress("eN", &_eN);

  if (IsFHC) { 
    duneobj->norm_s = (1e21/1.905e21);
  } else {
    duneobj->norm_s = (1e21/1.5e21);
  }
  duneobj->pot_s = (pot)/1e21;

  duneobj->nEvents = static_cast<int>(_data->GetEntries());

  duneobj->rw_yrec = new double[duneobj->nEvents];
  duneobj->rw_erec_lep = new double[duneobj->nEvents];
  duneobj->rw_erec_had = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_reco_q = new double[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 

  duneobj->rw_eRecoP = new double[duneobj->nEvents];
  duneobj->rw_eRecoPip = new double[duneobj->nEvents];
  duneobj->rw_eRecoPim = new double[duneobj->nEvents];
  duneobj->rw_eRecoPi0 = new double[duneobj->nEvents];
  duneobj->rw_eRecoN = new double[duneobj->nEvents];

  duneobj->rw_LepE = new double[duneobj->nEvents];
  duneobj->rw_eP = new double[duneobj->nEvents];
  duneobj->rw_ePip = new double[duneobj->nEvents];
  duneobj->rw_ePim = new double[duneobj->nEvents];
  duneobj->rw_ePi0 = new double[duneobj->nEvents];
  duneobj->rw_eN = new double[duneobj->nEvents];

  duneobj->nupdg = new int[duneobj->nEvents];
  duneobj->nupdgUnosc = new int[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);

  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    _data->GetEntry(i);

    duneobj->nupdg[i] = sample_nupdg[iSample];
    duneobj->nupdgUnosc[i] = sample_nupdgunosc[iSample];
    
    duneobj->rw_erec[i] = _erec;
    duneobj->rw_erec_shifted[i] = _erec;
    duneobj->rw_erec_lep[i] = _erec_lep;
    duneobj->rw_erec_had[i] = (_erec - _erec_lep);
    duneobj->rw_yrec[i] = ((_erec - _erec_lep)/_erec);
    duneobj->rw_etru[i] = _ev; // in GeV
    duneobj->rw_theta[i] = _LepNuAngle;
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_reco_q[i] = _reco_q;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;
    
    duneobj->rw_eRecoP[i] = _eRecoP; 
    duneobj->rw_eRecoPip[i] = _eRecoPip; 
    duneobj->rw_eRecoPim[i] = _eRecoPim; 
    duneobj->rw_eRecoPi0[i] = _eRecoPi0; 
    duneobj->rw_eRecoN[i] = _eRecoN; 
    
    duneobj->rw_LepE[i] = _LepE; 
    duneobj->rw_eP[i] = _eP; 
    duneobj->rw_ePip[i] = _ePip; 
    duneobj->rw_ePim[i] = _ePim; 
    duneobj->rw_ePi0[i] = _ePi0; 
    duneobj->rw_eN[i] = _eN; 
    
    //Assume everything is on Argon for now....
    duneobj->Target[i] = 40;

    int M3Mode = Modes->GetModeFromGenerator(std::abs(_mode));
    if (!_isCC) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    duneobj->mode[i] = M3Mode;
    
    duneobj->flux_w[i] = 1.0;
  }
  
  //_sampleFile->Close();
  _data->Reset();
  delete _data;
  return duneobj->nEvents;
}


const double* SampleHandlerBeamND::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunendmcSamples[iSample].rw_etru[iEvent]);
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &(dunendmcSamples[iSample].rw_erec_shifted[iEvent]);
    break;
  case kyRec:
    KinematicValue = &(dunendmcSamples[iSample].rw_yrec[iEvent]);
    break;
  case kOscChannel:
    KinematicValue = &(MCSamples[iSample].ChannelIndex);
    break;
  case kMode:
    KinematicValue = &(dunendmcSamples[iSample].mode[iEvent]);
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

const double* SampleHandlerBeamND::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* SampleHandlerBeamND::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double SampleHandlerBeamND::ReturnKinematicParameter(int KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iSample, iEvent);
}

double SampleHandlerBeamND::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iSample, iEvent);
}

void SampleHandlerBeamND::SetupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
    fdobj->nupdgUnosc[iEvent] = &(duneobj->nupdgUnosc[iEvent]);
    fdobj->nupdg[iEvent] = &(duneobj->nupdg[iEvent]);
  }
}


std::vector<double> SampleHandlerBeamND::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  std::vector<double> ReturnVec;

  switch (KinPar) {

  case kIsFHC:
    ReturnVec.resize(3);
    ReturnVec[0] = -0.5;
    ReturnVec[1] = 0.5;
    ReturnVec[2] = 1.5;
    break;
    
  case kTrueNeutrinoEnergy:
    for (int i=0;i<20;i++) {
      ReturnVec.emplace_back(i);
    }
    ReturnVec.emplace_back(100.);
    ReturnVec.emplace_back(1000.);
    break;

  case kRecoNeutrinoEnergy:
    ReturnVec.resize(XBinEdges.size());
    for (unsigned int bin_i=0;bin_i<XBinEdges.size();bin_i++) {ReturnVec[bin_i] = XBinEdges[bin_i];}
    break;

  case kyRec:
    ReturnVec.resize(YBinEdges.size());
    for (unsigned int bin_i=0;bin_i<YBinEdges.size();bin_i++) {ReturnVec[bin_i] = YBinEdges[bin_i];}
    break;

  case kOscChannel:
    ReturnVec.resize(GetNsamples());
    for (int bin_i=0;bin_i<GetNsamples();bin_i++) {ReturnVec[bin_i] = bin_i;}
    break;

  case kMode:
    ReturnVec.resize(Modes->GetNModes());
    for (int bin_i=0;bin_i<Modes->GetNModes();bin_i++) {ReturnVec[bin_i] = bin_i;}
    break;

  default:
    MACH3LOG_ERROR("Unknown KinPar: {}",static_cast<int>(KinPar));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return ReturnVec;
}

// Set the covariance matrix for this class
void SampleHandlerBeamND::setNDCovMatrix() {
  if (NDCovMatrix == NULL) {
    std::cerr << "Could not find ND Detector covariance matrix you provided to setCovMatrix" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  int nXBinsNDCov = static_cast<int>(XBinEdges.size()-1);
  int nYBinsNDCov = static_cast<int>(YBinEdges.size()-1);
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
        double CV = SampleHandlerFD_data[yBin][xBin];
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

  int nXBinsNDCov = static_cast<int>(XBinEdges.size()-1);
  int nYBinsNDCov = static_cast<int>(YBinEdges.size()-1);

  int covSize = nXBinsNDCov*nYBinsNDCov;

  std::vector<double> FlatData;
  std::vector<double> FlatMCPred;

  //2D -> 1D 
  for (int xBin = 0; xBin < nXBinsNDCov; xBin++) 
  {
    for (int yBin = 0; yBin < nYBinsNDCov; yBin++) 
    {
        double MCPred = SampleHandlerFD_array[yBin][xBin];
        FlatMCPred.push_back(MCPred);

        double DataVal = SampleHandlerFD_data[yBin][xBin];
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
