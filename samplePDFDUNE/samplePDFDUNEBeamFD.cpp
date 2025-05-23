#include "samplePDFDUNEBeamFD.h"

samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

samplePDFDUNEBeamFD::~samplePDFDUNEBeamFD() {
}

void samplePDFDUNEBeamFD::Init() {
  dunemcSamples.resize(nSamples,dunemc_base());
  
  if (CheckNodeExists(SampleManager->raw(), "DUNESampleBools", "iselike" )) {
    iselike = SampleManager->raw()["DUNESampleBools"]["iselike"].as<bool>();
  } else{
    MACH3LOG_ERROR("Did not find DUNESampleBools:iselike in {}, please add this", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (CheckNodeExists(SampleManager->raw(), "DUNESampleBools", "isFHC" )) {
    isFHC = SampleManager->raw()["DUNESampleBools"]["isFHC"].as<double>();
  } else{
    MACH3LOG_ERROR("Did not find DUNESampleBools:isFHC in {}, please add this", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  if (CheckNodeExists(SampleManager->raw(), "POT")) {
    pot = SampleManager->raw()["POT"].as<double>();
  } else{
    MACH3LOG_ERROR("POT not defined in {}, please add this!", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
}

void samplePDFDUNEBeamFD::SetupSplines() {

  ///@todo move all of the spline setup into core
  if(XsecCov->GetNumParamsFromSampleName(SampleName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = std::unique_ptr<splineFDBase>(new splinesDUNE(XsecCov,Modes));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = nullptr;
  }
  
  return;
}


// === HH: Functional parameters ===
void samplePDFDUNEBeamFD::TotalEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  // Total energy scale uncertainties for anything but CC Numu, see:
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39

  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_had[iEvent];
}

void samplePDFDUNEBeamFD::TotalEScaleNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent) {
  // A special case for Not (CC Numu), where we also scale Erec by lepton energy
  // Since we reconstruct muon energy in a different way, see:
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_lep[iEvent];
}

void samplePDFDUNEBeamFD::TotalEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_had[iEvent] * dunemcSamples[iSample].rw_erec_had_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::TotalEScaleSqrtNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in TotalEScaleNotCCNumu
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_lep[iEvent] * dunemcSamples[iSample].rw_erec_lep_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::TotalEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  // Erec/sqrt(Erec) = sqrt(Erec)
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_had_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::TotalEScaleInvSqrtNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in TotalEScaleNotCCNumu
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_lep_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::HadEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_sum_ehad[iEvent];
}

void samplePDFDUNEBeamFD::HadEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_sum_ehad[iEvent] * dunemcSamples[iSample].rw_sum_ehad_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::HadEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_sum_ehad_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::MuEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  // HH TODO: Functionally this is the same as TotalEScaleNotCCNumu, not sure if this function is even needed
  TotalEScaleNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamFD::MuEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in MuEScale
  TotalEScaleSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamFD::MuEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in MuEScale
  TotalEScaleInvSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamFD::NEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_eRecoN[iEvent];
}

void samplePDFDUNEBeamFD::NEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_eRecoN[iEvent] * dunemcSamples[iSample].rw_eRecoN_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::NEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_eRecoN_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::EMEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_eRecoPi0[iEvent];
}

void samplePDFDUNEBeamFD::EMEScaleCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // Again this is the same as TotalEScaleNotCCNumu, not sure if this function is needed
  TotalEScaleNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamFD::EMEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_eRecoPi0[iEvent] * dunemcSamples[iSample].rw_eRecoPi0_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::EMEScaleSqrtCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in EMEScaleCCNue
  TotalEScaleSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamFD::EMEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_eRecoPi0_sqrt[iEvent];
}

void samplePDFDUNEBeamFD::EMEScaleInvSqrtCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in EMEScaleCCNue
  TotalEScaleInvSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamFD::HadRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True sum - reco sum
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunemcSamples[iSample].rw_eP[iEvent]
    + dunemcSamples[iSample].rw_ePip[iEvent]
    + dunemcSamples[iSample].rw_ePim[iEvent]
    - dunemcSamples[iSample].rw_sum_ehad[iEvent]);
}

void samplePDFDUNEBeamFD::MuRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True muon energy - reco muon energy
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunemcSamples[iSample].rw_LepE[iEvent] - dunemcSamples[iSample].rw_erec_lep[iEvent]);
}

void samplePDFDUNEBeamFD::NRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True neutron energy - reco neutron energy
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunemcSamples[iSample].rw_eN[iEvent] - dunemcSamples[iSample].rw_eRecoN[iEvent]);
}

void samplePDFDUNEBeamFD::EMRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True pi0 energy - reco pi0 energy
  dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunemcSamples[iSample].rw_ePi0[iEvent] - dunemcSamples[iSample].rw_eRecoPi0[iEvent]);
}

void samplePDFDUNEBeamFD::EMResCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // This is the same as MuRes, again not sure if this function is needed
  MuRes(par, iSample, iEvent);
}

void samplePDFDUNEBeamFD::RecoCVNNumu(const double * par, std::size_t iSample, std::size_t iEvent) {
  // CVN numu uncertainty
  dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent] += (*par);
}

void samplePDFDUNEBeamFD::RecoCVNNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // CVN nue uncertainty
  dunemcSamples[iSample].rw_cvnnue_shifted[iEvent] += (*par);
}

void samplePDFDUNEBeamFD::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering functional parameters");
  // This function manually populates the map of functional parameters
  // Maps the name of the functional parameter to the pointer of the function

  // This is the part where we manually enter things
  // A lambda function has to be used so we can refer to a non-static member function
  RegisterIndividualFuncPar("TotalEScaleFD",
                            kTotalEScale,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleNotCCNumuFD",
                            kTotalEScaleNotCCNumu,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleNotCCNumu(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleSqrtFD",
                            kTotalEScaleSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleSqrtNotCCNumuFD",
                            kTotalEScaleSqrtNotCCNumu,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleSqrtNotCCNumu(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleInvSqrtFD",
                            kTotalEScaleInvSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleInvSqrtNotCCNumuFD",
                            kTotalEScaleInvSqrtNotCCNumu,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleInvSqrtNotCCNumu(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadEScaleFD",
                            kHadEScale,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadEScaleSqrtFD",
                            kHadEScaleSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadEScaleInvSqrtFD",
                            kHadEScaleInvSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuEScaleFD",
                            kMuEScale,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuEScaleSqrtFD",
                            kMuEScaleSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuEScaleInvSqrtFD",
                            kMuEScaleInvSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NEScaleFD",
                            kNEScale,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NEScaleSqrtFD",
                            kNEScaleSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NEScaleInvSqrtFD",
                            kNEScaleInvSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleFD",
                            kEMEScale,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleCCNueFD",
                            kEMEScaleCCNue,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleCCNue(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleSqrtFD",
                            kEMEScaleSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleSqrtCCNueFD",
                            kEMEScaleSqrtCCNue,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleSqrtCCNue(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleInvSqrtFD",
                            kEMEScaleInvSqrt,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleInvSqrtCCNueFD",
                            kEMEScaleInvSqrtCCNue,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleInvSqrtCCNue(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadResFD",
                            kHadRes,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuResFD",
                            kMuRes,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NResFD",
                            kNRes,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMResFD",
                            kEMRes,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMResCCNueFD",
                            kEMResCCNue,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMResCCNue(par, iSample, iEvent); });

  RegisterIndividualFuncPar("RecoCVNNumuFD",
                            kRecoCVNNumu,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->RecoCVNNumu(par, iSample, iEvent); });

  RegisterIndividualFuncPar("RecoCVNNueFD",
                            kRecoCVNNue,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->RecoCVNNue(par, iSample, iEvent); });

  MACH3LOG_INFO("Finished registering functional parameters");
}

// HH: Reset the shifted values to the original values
void samplePDFDUNEBeamFD::resetShifts(int iSample, int iEvent) {
  dunemcSamples[iSample].rw_erec_shifted[iEvent] = dunemcSamples[iSample].rw_erec[iEvent];
  dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnumu[iEvent];
  dunemcSamples[iSample].rw_cvnnue_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnue[iEvent];
}
// =================================

void samplePDFDUNEBeamFD::SetupWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunemcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunemcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunemcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}


int samplePDFDUNEBeamFD::setupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files[iSample]);
  
  TFile* _sampleFile = TFile::Open(mc_files[iSample].c_str(), "READ");
  TTree* _data = _sampleFile->Get<TTree>("caf");
  
  if(_data){
    MACH3LOG_INFO("Found \"caf\" tree in {}", mc_files[iSample]);
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find \"caf\" tree in {}", mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //Reco Variables
  double _erec;
  double _erec_nue;
  double _erec_had;
  double _erec_had_nue;
  double _erec_lep;
  double _erec_lep_nue;

  double _eRecoP;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoN;

  double _cvnnumu;
  double _cvnnue;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;

  //Truth Variables
  int _mode;  
  double _ev;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _BeRPA_cvwgt;
  int _isCC;
  int _nuPDGunosc;
  int _nuPDG;
  
  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
  _data->SetBranchStatus("Ev_reco_numu", 1);
  _data->SetBranchAddress("Ev_reco_numu", &_erec);
  _data->SetBranchStatus("Ev_reco_nue", 1);
  _data->SetBranchAddress("Ev_reco_nue", &_erec_nue);
  _data->SetBranchStatus("RecoHadEnNumu", 1);
  _data->SetBranchAddress("RecoHadEnNumu", &_erec_had);
  _data->SetBranchStatus("RecoHadEnNue", 1);
  _data->SetBranchAddress("RecoHadEnNue", &_erec_had_nue);
  _data->SetBranchStatus("RecoLepEnNumu", 1);
  _data->SetBranchAddress("RecoLepEnNumu", &_erec_lep);
  _data->SetBranchStatus("RecoLepEnNue", 1);
  _data->SetBranchAddress("RecoLepEnNue", &_erec_lep_nue);

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

  _data->SetBranchStatus("mode",1);
  _data->SetBranchAddress("mode",&_mode);
  _data->SetBranchStatus("cvnnumu",1);
  _data->SetBranchAddress("cvnnumu", &_cvnnumu);
  _data->SetBranchStatus("cvnnue",1);
  _data->SetBranchAddress("cvnnue", &_cvnnue);
  _data->SetBranchStatus("isCC", 1);
  _data->SetBranchAddress("isCC", &_isCC);
  _data->SetBranchStatus("nuPDGunosc", 1);
  _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
  _data->SetBranchStatus("nuPDG", 1);
  _data->SetBranchAddress("nuPDG", &_nuPDG);
  _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
  _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);
  _data->SetBranchStatus("vtx_x", 1);
  _data->SetBranchAddress("vtx_x", &_vtx_x);
  _data->SetBranchStatus("vtx_y", 1);
  _data->SetBranchAddress("vtx_y", &_vtx_y);
  _data->SetBranchStatus("vtx_z", 1);
  _data->SetBranchAddress("vtx_z", &_vtx_z);  

  TH1D* norm = _sampleFile->Get<TH1D>("norm");
  if(!norm){
    MACH3LOG_ERROR("Add a norm KEY to the root file using MakeNormHists.cxx");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // now fill the actual variables
  duneobj->norm_s = norm->GetBinContent(1);
  duneobj->pot_s = pot/norm->GetBinContent(2);

  duneobj->nEvents = static_cast<int>(_data->GetEntries());

  // allocate memory for dunemc variables
  duneobj->rw_cvnnumu = new double[duneobj->nEvents];
  duneobj->rw_cvnnue = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu_shifted = new double[duneobj->nEvents];
  duneobj->rw_cvnnue_shifted = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_erec_had = new double[duneobj->nEvents];
  duneobj->rw_erec_lep = new double[duneobj->nEvents];

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

  duneobj->rw_erec_had_sqrt = new double[duneobj->nEvents];
  duneobj->rw_erec_lep_sqrt = new double[duneobj->nEvents];
  duneobj->rw_eRecoN_sqrt = new double[duneobj->nEvents];
  duneobj->rw_eRecoPi0_sqrt = new double[duneobj->nEvents];

  duneobj->rw_sum_ehad = new double[duneobj->nEvents];
  duneobj->rw_sum_ehad_sqrt = new double[duneobj->nEvents];

  duneobj->rw_trueccnue = new double[duneobj->nEvents];
  duneobj->rw_trueccnumu = new double[duneobj->nEvents];


  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];

  duneobj->nupdgUnosc = new int[duneobj->nEvents];
  duneobj->nupdg = new int[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);
  
  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    _data->GetEntry(i);

    duneobj->nupdg[i] = sample_nupdg[iSample];
    duneobj->nupdgUnosc[i] = sample_nupdgunosc[iSample];    
    
    duneobj->rw_cvnnumu[i] = (_cvnnumu);
    duneobj->rw_cvnnue[i] = (_cvnnue);
    duneobj->rw_cvnnumu_shifted[i] = (_cvnnumu); 
    duneobj->rw_cvnnue_shifted[i] = (_cvnnue);
    if (iselike) {
      duneobj->rw_erec[i] = (_erec_nue);
      duneobj->rw_erec_shifted[i] = (_erec_nue); 
      duneobj->rw_erec_had[i] = (_erec_had_nue);
      duneobj->rw_erec_lep[i] = (_erec_lep_nue);
    } else {
      duneobj->rw_erec[i] = (_erec); 
      duneobj->rw_erec_shifted[i] = (_erec); 
      duneobj->rw_erec_had[i] = (_erec_had); 
      duneobj->rw_erec_lep[i] = (_erec_lep); 
    }
    
    duneobj->rw_eRecoP[i] = (_eRecoP); 
    duneobj->rw_eRecoPip[i] = (_eRecoPip); 
    duneobj->rw_eRecoPim[i] = (_eRecoPim); 
    duneobj->rw_eRecoPi0[i] = (_eRecoPi0); 
    duneobj->rw_eRecoN[i] = (_eRecoN); 
    
    duneobj->rw_LepE[i] = (_LepE); 
    duneobj->rw_eP[i] = (_eP); 
    duneobj->rw_ePip[i] = (_ePip); 
    duneobj->rw_ePim[i] = (_ePim); 
    duneobj->rw_ePi0[i] = (_ePi0); 
    duneobj->rw_eN[i] = (_eN);

	duneobj->rw_erec_had_sqrt[i] = sqrt(duneobj->rw_erec_had[i]);
    duneobj->rw_erec_lep_sqrt[i] = sqrt(duneobj->rw_erec_lep[i]);
    duneobj->rw_eRecoN_sqrt[i] = sqrt(duneobj->rw_eRecoN[i]);
    duneobj->rw_eRecoPi0_sqrt[i] = sqrt(duneobj->rw_eRecoPi0[i]);

    duneobj->rw_sum_ehad[i] = duneobj->rw_eRecoP[i] + duneobj->rw_eRecoPip[i] + duneobj->rw_eRecoPim[i];
    duneobj->rw_sum_ehad_sqrt[i] = sqrt(duneobj->rw_sum_ehad[i]);

    duneobj->rw_etru[i] = (_ev);
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = (_BeRPA_cvwgt);
    duneobj->rw_vtx_x[i] = (_vtx_x);
    duneobj->rw_vtx_y[i] = (_vtx_y);
    duneobj->rw_vtx_z[i] = (_vtx_z);

    duneobj->rw_trueccnumu[i] = static_cast<double>(duneobj->rw_isCC[i]==1 && abs(duneobj->rw_nuPDG[i])==14);
    duneobj->rw_trueccnue[i] = static_cast<double>(duneobj->rw_isCC[i]==1 && abs(duneobj->rw_nuPDG[i])==12);
    
    //Assume everything is on Argon40 for now....
    duneobj->Target[i] = kTarget_Ar;
    
    int M3Mode = Modes->GetModeFromGenerator(std::abs(_mode));
    if (!_isCC) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    duneobj->mode[i] = M3Mode;
    
    duneobj->flux_w[i] = 1.0;
  }
  
  _sampleFile->Close();
  return duneobj->nEvents;
}

const double* samplePDFDUNEBeamFD::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue = nullptr;

  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &dunemcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &dunemcSamples[iSample].rw_erec_shifted[iEvent];
    break;
  case kTrueXPos:
    KinematicValue = &dunemcSamples[iSample].rw_vtx_x[iEvent];
    break;
  case kTrueYPos:
    KinematicValue = &dunemcSamples[iSample].rw_vtx_y[iEvent];
    break;
  case kTrueZPos:
    KinematicValue = &dunemcSamples[iSample].rw_vtx_z[iEvent];
    break;
  case kCVNNumu:
    KinematicValue = &dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
    break;
  case kCVNNue:
    KinematicValue = &dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
    break;
  case kM3Mode:
    KinematicValue = &(dunemcSamples[iSample].mode[iEvent]);
    break;
  case kOscChannel:
    KinematicValue = &(MCSamples[iSample].ChannelIndex);
    break;
  case kIsFHC:
    KinematicValue = &(isFHC);
    break;
  case kTrueCCnue: 
	KinematicValue = &(dunemcSamples[iSample].rw_trueccnue[iEvent]);
 	break;
  case kTrueCCnumu: 
	KinematicValue = &(dunemcSamples[iSample].rw_trueccnumu[iEvent]);
 	break;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* samplePDFDUNEBeamFD::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter)); 
  return GetPointerToKinematicParameter(KinPar, iSample, iEvent);
}

const double* samplePDFDUNEBeamFD::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar, iSample, iEvent);
}

double samplePDFDUNEBeamFD::ReturnKinematicParameter(int KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iSample, iEvent);
}

double samplePDFDUNEBeamFD::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
 return *GetPointerToKinematicParameter(KinematicParameter, iSample, iEvent);
}

void samplePDFDUNEBeamFD::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);  
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
    fdobj->nupdg[iEvent] = &(duneobj->nupdg[iEvent]);
    fdobj->nupdgUnosc[iEvent] = &(duneobj->nupdgUnosc[iEvent]);
  }
  
}
 
std::vector<double> samplePDFDUNEBeamFD::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinematicParameter = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  std::vector<double> ReturnVec;
  
  switch(KinematicParameter){
  case kIsFHC:
    ReturnVec.resize(3);
    ReturnVec[0] = -0.5;
    ReturnVec[1] = 0.5;
    ReturnVec[2] = 1.5;
    break;
    
  case kTrueNeutrinoEnergy:
  case kRecoNeutrinoEnergy:
    ReturnVec.resize(XBinEdges.size());
    for (unsigned int bin_i=0;bin_i<XBinEdges.size();bin_i++) {ReturnVec[bin_i] = XBinEdges[bin_i];}
    break;

  case kOscChannel:
    ReturnVec.resize(GetNsamples());
    for (int bin_i=0;bin_i<GetNsamples();bin_i++) {ReturnVec[bin_i] = bin_i;}
    break;

  case kM3Mode:
    ReturnVec.resize(Modes->GetNModes());
    for (int bin_i=0;bin_i<Modes->GetNModes();bin_i++) {ReturnVec[bin_i] = bin_i;}
    break;

  case kTrueXPos:
  case kTrueYPos:
  case kTrueZPos:
  case kTrueCCnue:
  case kTrueCCnumu:
  case kCVNNue:
  case kCVNNumu:
    ReturnVec.resize(2);
    ReturnVec[0] = 1e-8;
    ReturnVec[1] = 1e8;
    break;

  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", static_cast<int>(KinematicParameter));
    throw MaCh3Exception(__FILE__, __LINE__);

  }      
  
  return ReturnVec;
}
