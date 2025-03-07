#include "samplePDFDUNEBeamND.h"

samplePDFDUNEBeamND::samplePDFDUNEBeamND(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_=nullptr) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {
  OscCov = nullptr;
  
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

samplePDFDUNEBeamND::~samplePDFDUNEBeamND() {
}

void samplePDFDUNEBeamND::Init() {
  dunendmcSamples.resize(nSamples,dunemc_base());
  
  IsFHC = SampleManager->raw()["SampleBools"]["isFHC"].as<double>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();
  pot = SampleManager->raw()["POT"].as<double>();

  tot_escale_nd_pos = -999;
  tot_escale_sqrt_nd_pos = -999;
  tot_escale_invsqrt_nd_pos = -999;
  had_escale_nd_pos = -999;
  had_escale_sqrt_nd_pos = -999;
  had_escale_invsqrt_nd_pos = -999;
  mu_escale_nd_pos = -999;
  mu_escale_sqrt_nd_pos = -999;
  mu_escale_invsqrt_nd_pos = -999;
  n_escale_nd_pos = -999;
  n_escale_sqrt_nd_pos = -999;
  n_escale_invsqrt_nd_pos = -999;
  em_escale_nd_pos = -999;
  em_escale_sqrt_nd_pos = -999;
  em_escale_invsqrt_nd_pos = -999;
  had_res_nd_pos = -999;
  mu_res_nd_pos = -999;
  n_res_nd_pos = -999;
  em_res_nd_pos = -999;

  /*
  nNDDetectorSystPointers = funcParsIndex.size();
  NDDetectorSystPointers = std::vector<const double*>(nNDDetectorSystPointers);

  int func_it = 0;
  for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++func_it) {
    std::string name = funcParsNames.at(func_it);
    
    if (name == "TotalEScaleND") {
      tot_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_nd_pos);
    }
    else if (name == "TotalEScaleSqrtND") {
      tot_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_sqrt_nd_pos);
    }
    else if (name == "TotalEScaleInvSqrtND") {
      tot_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_invsqrt_nd_pos);
    }
    else if (name == "HadEScaleND") {
      had_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_nd_pos);
    }
    else if (name == "HadEScaleSqrtND") {
      had_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_sqrt_nd_pos);
    }
    else if (name == "HadEScaleInvSqrtND") {
      had_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_invsqrt_nd_pos);
    }
    else if (name == "MuEScaleND") {
      mu_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_nd_pos);
    }
    else if (name == "MuEScaleSqrtND") {
      mu_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_sqrt_nd_pos);
    }
    else if (name == "MuEScaleInvSqrtND") {
      mu_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_invsqrt_nd_pos);
    }
    else if (name == "NEScaleND") {
      n_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_nd_pos);
    }
    else if (name == "NEScaleSqrtND") {
      n_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_sqrt_nd_pos);
    }
    else if (name == "NEScaleInvSqrtND") {
      n_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_invsqrt_nd_pos);
    }
    else if (name == "EMEScaleND") {
      em_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_nd_pos);
    }
    else if (name == "EMEScaleSqrtND") {
      em_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_sqrt_nd_pos);
    }
    else if (name == "EMEScaleInvSqrtND") {
      em_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_invsqrt_nd_pos);
    }
    else if (name == "HadResND") {
      had_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_res_nd_pos);
    }
    else if (name == "MuResND") {
      mu_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_res_nd_pos);
    }
    else if (name == "NResND") {
      n_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_res_nd_pos);
    }
    else if (name == "EMResND") {
      em_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_res_nd_pos);
    }
    else { 
      MACH3LOG_ERROR("Found a functional parameter which wasn't specified in the xml | samplePDFDUNEBeamND: {}",name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  */
  
  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void samplePDFDUNEBeamND::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    SplineHandler = std::unique_ptr<splineFDBase>(new splinesDUNE(XsecCov,Modes));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

void samplePDFDUNEBeamND::SetupWeightPointers() {
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

int samplePDFDUNEBeamND::setupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files.at(iSample));
  
  _sampleFile = TFile::Open(mc_files.at(iSample).c_str(), "READ");
  _data = _sampleFile->Get<TTree>("caf");

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
    duneobj->norm_s = (1e21/1.5e21);
  } else {
    duneobj->norm_s = (1e21/1.905e21);
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
  
  _sampleFile->Close();
  return duneobj->nEvents;
}


const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
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

const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamND::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicVariable, iSample, iEvent);
}

double samplePDFDUNEBeamND::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iSample, iEvent);
}

void samplePDFDUNEBeamND::setupFDMC(int iSample) {
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

void samplePDFDUNEBeamND::applyShifts(int iSample, int iEvent) {
  // reset erec back to original value
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] = dunendmcSamples[iSample].rw_erec[iEvent];

  /*
  //Calculate values needed
  double sqrtErecHad =  sqrt(dunendmcSamples[iSample].rw_erec_had[iEvent]);
  double sqrtErecLep =  sqrt(dunendmcSamples[iSample].rw_erec_lep[iEvent]);
  double sqrteRecoPi0 = sqrt(dunendmcSamples[iSample].rw_eRecoPi0[iEvent]);
  double sqrteRecoN = sqrt(dunendmcSamples[iSample].rw_eRecoN[iEvent]);
  double sumEhad = dunendmcSamples[iSample].rw_eRecoP[iEvent] + dunendmcSamples[iSample].rw_eRecoPip[iEvent] + dunendmcSamples[iSample].rw_eRecoPim[iEvent];
  double sqrtSumEhad = sqrt(sumEhad);

  double invSqrtErecHad =  1/(sqrtErecHad+0.1);
  double invSqrtErecLep =  1/(sqrtErecLep+0.1);
  double invSqrteRecoPi0 =  1/(sqrteRecoPi0+0.1);
  double invSqrteRecoN =  1/(sqrteRecoN+0.1);
  double invSqrtSumEhad =  1/(sqrtSumEhad+0.1);

  bool CCnumu {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==14 && dunendmcSamples[iSample].nupdgUnosc[iEvent]==2};
  bool CCnue {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==12 && dunendmcSamples[iSample].nupdgUnosc[iEvent]==1};
  bool NotCCnumu {!(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==14) && dunendmcSamples[iSample].nupdgUnosc[iEvent]==2};


  TotalEScaleND(NDDetectorSystPointers[0], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_had[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], NotCCnumu);

  TotalEScaleSqrtND(NDDetectorSystPointers[1], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_had[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], sqrtErecHad, sqrtErecLep, NotCCnumu);

  TotalEScaleInvSqrtND(NDDetectorSystPointers[2], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_had[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecHad, invSqrtErecLep, NotCCnumu);

  HadEScaleND(NDDetectorSystPointers[3], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], sumEhad);

  HadEScaleSqrtND(NDDetectorSystPointers[4], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, sqrtSumEhad);

  HadEScaleInvSqrtND(NDDetectorSystPointers[5], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, invSqrtSumEhad);

  MuEScaleND(NDDetectorSystPointers[6], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], CCnumu);

  MuEScaleSqrtND(NDDetectorSystPointers[7], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, CCnumu);

  MuEScaleInvSqrtND(NDDetectorSystPointers[8], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, CCnumu);

  NEScaleND(NDDetectorSystPointers[9], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent]);

  NEScaleSqrtND(NDDetectorSystPointers[10], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent], sqrteRecoN);

  NEScaleInvSqrtND(NDDetectorSystPointers[11], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent], invSqrteRecoN);

  EMEScaleND(NDDetectorSystPointers[12], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], CCnue);

  EMEScaleSqrtND(NDDetectorSystPointers[13], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, sqrteRecoPi0, CCnue);

  EMEScaleInvSqrtND(NDDetectorSystPointers[14], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, invSqrteRecoPi0, CCnue);

  HadResND(NDDetectorSystPointers[15], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoP[iEvent], dunendmcSamples[iSample].rw_eRecoPip[iEvent], dunendmcSamples[iSample].rw_eRecoPim[iEvent], dunendmcSamples[iSample].rw_eP[iEvent], dunendmcSamples[iSample].rw_ePip[iEvent], dunendmcSamples[iSample].rw_ePim[iEvent]);

  MuResND(NDDetectorSystPointers[16], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], dunendmcSamples[iSample].rw_LepE[iEvent], CCnumu);

  NResND(NDDetectorSystPointers[17], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent], dunendmcSamples[iSample].rw_eN[iEvent]);

  EMResND(NDDetectorSystPointers[18], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_ePi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], dunendmcSamples[iSample].rw_LepE[iEvent], CCnue);
  */
}

std::vector<double> samplePDFDUNEBeamND::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
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
