#include "SampleHandlerBeamFD.h"

SampleHandlerBeamFD::SampleHandlerBeamFD(std::string mc_version_, ParameterHandlerGeneric* ParHandler_) : SampleHandlerFD(mc_version_, ParHandler_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerBeamFD::~SampleHandlerBeamFD() {
}

void SampleHandlerBeamFD::Init() {
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

void SampleHandlerBeamFD::SetupSplines() {

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


// === HH: Functional parameters ===
void SampleHandlerBeamFD::TotalEScale(const double * par, std::size_t iEvent) {
  // Total energy scale uncertainties for anything but CC Numu, see:
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39

  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_erec_had;
}

void SampleHandlerBeamFD::TotalEScaleNotCCNumu(const double * par, std::size_t iEvent) {
  // A special case for Not (CC Numu), where we also scale Erec by lepton energy
  // Since we reconstruct muon energy in a different way, see:
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_erec_lep;
}

void SampleHandlerBeamFD::TotalEScaleSqrt(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_erec_had * dunemcSamples[iEvent].rw_erec_had_sqrt;
}

void SampleHandlerBeamFD::TotalEScaleSqrtNotCCNumu(const double * par, std::size_t iEvent) {
  // See comments in TotalEScaleNotCCNumu
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_erec_lep * dunemcSamples[iEvent].rw_erec_lep_sqrt;
}

void SampleHandlerBeamFD::TotalEScaleInvSqrt(const double * par, std::size_t iEvent) {
  // Erec/sqrt(Erec) = sqrt(Erec)
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_erec_had_sqrt;
}

void SampleHandlerBeamFD::TotalEScaleInvSqrtNotCCNumu(const double * par, std::size_t iEvent) {
  // See comments in TotalEScaleNotCCNumu
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_erec_lep_sqrt;
}

void SampleHandlerBeamFD::HadEScale(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_sum_ehad;
}

void SampleHandlerBeamFD::HadEScaleSqrt(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_sum_ehad * dunemcSamples[iEvent].rw_sum_ehad_sqrt;
}

void SampleHandlerBeamFD::HadEScaleInvSqrt(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_sum_ehad_sqrt;
}

void SampleHandlerBeamFD::MuEScale(const double * par, std::size_t iEvent) {
  // HH TODO: Functionally this is the same as TotalEScaleNotCCNumu, not sure if this function is even needed
  TotalEScaleNotCCNumu(par, iEvent);
}

void SampleHandlerBeamFD::MuEScaleSqrt(const double * par, std::size_t iEvent) {
  // See comments in MuEScale
  TotalEScaleSqrtNotCCNumu(par, iEvent);
}

void SampleHandlerBeamFD::MuEScaleInvSqrt(const double * par, std::size_t iEvent) {
  // See comments in MuEScale
  TotalEScaleInvSqrtNotCCNumu(par, iEvent);
}

void SampleHandlerBeamFD::NEScale(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_eRecoN;
}

void SampleHandlerBeamFD::NEScaleSqrt(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_eRecoN * dunemcSamples[iEvent].rw_eRecoN_sqrt;
}

void SampleHandlerBeamFD::NEScaleInvSqrt(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_eRecoN_sqrt;
}

void SampleHandlerBeamFD::EMEScale(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_eRecoPi0;
}

void SampleHandlerBeamFD::EMEScaleCCNue(const double * par, std::size_t iEvent) {
  // Again this is the same as TotalEScaleNotCCNumu, not sure if this function is needed
  TotalEScaleNotCCNumu(par, iEvent);
}

void SampleHandlerBeamFD::EMEScaleSqrt(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_eRecoPi0 * dunemcSamples[iEvent].rw_eRecoPi0_sqrt;
}

void SampleHandlerBeamFD::EMEScaleSqrtCCNue(const double * par, std::size_t iEvent) {
  // See comments in EMEScaleCCNue
  TotalEScaleSqrtNotCCNumu(par, iEvent);
}

void SampleHandlerBeamFD::EMEScaleInvSqrt(const double * par, std::size_t iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * dunemcSamples[iEvent].rw_eRecoPi0_sqrt;
}

void SampleHandlerBeamFD::EMEScaleInvSqrtCCNue(const double * par, std::size_t iEvent) {
  // See comments in EMEScaleCCNue
  TotalEScaleInvSqrtNotCCNumu(par, iEvent);
}

void SampleHandlerBeamFD::HadRes(const double * par, std::size_t iEvent) {
  // True sum - reco sum
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * (dunemcSamples[iEvent].rw_eP
    + dunemcSamples[iEvent].rw_ePip
    + dunemcSamples[iEvent].rw_ePim
    - dunemcSamples[iEvent].rw_sum_ehad);
}

void SampleHandlerBeamFD::MuRes(const double * par, std::size_t iEvent) {
  // True muon energy - reco muon energy
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * (dunemcSamples[iEvent].rw_LepE - dunemcSamples[iEvent].rw_erec_lep);
}

void SampleHandlerBeamFD::NRes(const double * par, std::size_t iEvent) {
  // True neutron energy - reco neutron energy
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * (dunemcSamples[iEvent].rw_eN - dunemcSamples[iEvent].rw_eRecoN);
}

void SampleHandlerBeamFD::EMRes(const double * par, std::size_t iEvent) {
  // True pi0 energy - reco pi0 energy
  dunemcSamples[iEvent].rw_erec_shifted += (*par) * (dunemcSamples[iEvent].rw_ePi0 - dunemcSamples[iEvent].rw_eRecoPi0);
}

void SampleHandlerBeamFD::EMResCCNue(const double * par, std::size_t iEvent) {
  // This is the same as MuRes, again not sure if this function is needed
  MuRes(par, iEvent);
}

void SampleHandlerBeamFD::RecoCVNNumu(const double * par, std::size_t iEvent) {
  // CVN numu uncertainty
  dunemcSamples[iEvent].rw_cvnnumu_shifted += (*par);
}

void SampleHandlerBeamFD::RecoCVNNue(const double * par, std::size_t iEvent) {
  // CVN nue uncertainty
  dunemcSamples[iEvent].rw_cvnnue_shifted += (*par);
}

void SampleHandlerBeamFD::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering functional parameters");
  // This function manually populates the map of functional parameters
  // Maps the name of the functional parameter to the pointer of the function

  // This is the part where we manually enter things
  // A lambda function has to be used so we can refer to a non-static member function
  RegisterIndividualFunctionalParameter("TotalEScaleFD",
                            kTotalEScale,
                            [this](const double * par, std::size_t iEvent) { this->TotalEScale(par, iEvent); });

  RegisterIndividualFunctionalParameter("TotalEScaleNotCCNumuFD",
                            kTotalEScaleNotCCNumu,
                            [this](const double * par, std::size_t iEvent) { this->TotalEScaleNotCCNumu(par, iEvent); });

  RegisterIndividualFunctionalParameter("TotalEScaleSqrtFD",
                            kTotalEScaleSqrt,
                            [this](const double * par, std::size_t iEvent) { this->TotalEScaleSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("TotalEScaleSqrtNotCCNumuFD",
                            kTotalEScaleSqrtNotCCNumu,
                            [this](const double * par, std::size_t iEvent) { this->TotalEScaleSqrtNotCCNumu(par, iEvent); });

  RegisterIndividualFunctionalParameter("TotalEScaleInvSqrtFD",
                            kTotalEScaleInvSqrt,
                            [this](const double * par, std::size_t iEvent) { this->TotalEScaleInvSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("TotalEScaleInvSqrtNotCCNumuFD",
                            kTotalEScaleInvSqrtNotCCNumu,
                            [this](const double * par, std::size_t iEvent) { this->TotalEScaleInvSqrtNotCCNumu(par, iEvent); });

  RegisterIndividualFunctionalParameter("HadEScaleFD",
                            kHadEScale,
                            [this](const double * par, std::size_t iEvent) { this->HadEScale(par, iEvent); });

  RegisterIndividualFunctionalParameter("HadEScaleSqrtFD",
                            kHadEScaleSqrt,
                            [this](const double * par, std::size_t iEvent) { this->HadEScaleSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("HadEScaleInvSqrtFD",
                            kHadEScaleInvSqrt,
                            [this](const double * par, std::size_t iEvent) { this->HadEScaleInvSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("MuEScaleFD",
                            kMuEScale,
                            [this](const double * par, std::size_t iEvent) { this->MuEScale(par, iEvent); });

  RegisterIndividualFunctionalParameter("MuEScaleSqrtFD",
                            kMuEScaleSqrt,
                            [this](const double * par, std::size_t iEvent) { this->MuEScaleSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("MuEScaleInvSqrtFD",
                            kMuEScaleInvSqrt,
                            [this](const double * par, std::size_t iEvent) { this->MuEScaleInvSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("NEScaleFD",
                            kNEScale,
                            [this](const double * par, std::size_t iEvent) { this->NEScale(par, iEvent); });

  RegisterIndividualFunctionalParameter("NEScaleSqrtFD",
                            kNEScaleSqrt,
                            [this](const double * par, std::size_t iEvent) { this->NEScaleSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("NEScaleInvSqrtFD",
                            kNEScaleInvSqrt,
                            [this](const double * par, std::size_t iEvent) { this->NEScaleInvSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMEScaleFD",
                            kEMEScale,
                            [this](const double * par, std::size_t iEvent) { this->EMEScale(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMEScaleCCNueFD",
                            kEMEScaleCCNue,
                            [this](const double * par, std::size_t iEvent) { this->EMEScaleCCNue(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMEScaleSqrtFD",
                            kEMEScaleSqrt,
                            [this](const double * par, std::size_t iEvent) { this->EMEScaleSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMEScaleSqrtCCNueFD",
                            kEMEScaleSqrtCCNue,
                            [this](const double * par, std::size_t iEvent) { this->EMEScaleSqrtCCNue(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMEScaleInvSqrtFD",
                            kEMEScaleInvSqrt,
                            [this](const double * par, std::size_t iEvent) { this->EMEScaleInvSqrt(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMEScaleInvSqrtCCNueFD",
                            kEMEScaleInvSqrtCCNue,
                            [this](const double * par, std::size_t iEvent) { this->EMEScaleInvSqrtCCNue(par, iEvent); });

  RegisterIndividualFunctionalParameter("HadResFD",
                            kHadRes,
                            [this](const double * par, std::size_t iEvent) { this->HadRes(par, iEvent); });

  RegisterIndividualFunctionalParameter("MuResFD",
                            kMuRes,
                            [this](const double * par, std::size_t iEvent) { this->MuRes(par, iEvent); });

  RegisterIndividualFunctionalParameter("NResFD",
                            kNRes,
                            [this](const double * par, std::size_t iEvent) { this->NRes(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMResFD",
                            kEMRes,
                            [this](const double * par, std::size_t iEvent) { this->EMRes(par, iEvent); });

  RegisterIndividualFunctionalParameter("EMResCCNueFD",
                            kEMResCCNue,
                            [this](const double * par, std::size_t iEvent) { this->EMResCCNue(par, iEvent); });

  RegisterIndividualFunctionalParameter("RecoCVNNumuFD",
                            kRecoCVNNumu,
                            [this](const double * par, std::size_t iEvent) { this->RecoCVNNumu(par, iEvent); });

  RegisterIndividualFunctionalParameter("RecoCVNNueFD",
                            kRecoCVNNue,
                            [this](const double * par, std::size_t iEvent) { this->RecoCVNNue(par, iEvent); });

  MACH3LOG_INFO("Finished registering functional parameters");
}

// HH: Reset the shifted values to the original values
void SampleHandlerBeamFD::resetShifts(int iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted = dunemcSamples[iEvent].rw_erec;
  dunemcSamples[iEvent].rw_cvnnumu_shifted = dunemcSamples[iEvent].rw_cvnnumu;
  dunemcSamples[iEvent].rw_cvnnue_shifted = dunemcSamples[iEvent].rw_cvnnue;
}
// =================================

void SampleHandlerBeamFD::SetupWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    // MCSamples[i].ntotal_weight_pointers[j] = 6;
    // MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].pot_s));
    MCSamples[i].total_weight_pointers.push_back( &(dunemcSamples[i].norm_s));
    MCSamples[i].total_weight_pointers.push_back( MCSamples[i].osc_w_pointer);
    MCSamples[i].total_weight_pointers.push_back( &(dunemcSamples[i].rw_berpaacvwgt));
    MCSamples[i].total_weight_pointers.push_back( &(dunemcSamples[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back( &(MCSamples[i].xsec_w));
  }
}


int SampleHandlerBeamFD::SetupExperimentMC() {

  // dunemc_base *duneobj = &(dunemcSamples[iSample]);
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
  TChain* _data = new TChain("caf");
  for (size_t iSample=0;iSample<mc_files.size();iSample++) {
    MACH3LOG_INFO("Adding file to TChain: {}", mc_files[iSample]);
    TFile* _sampleFile = TFile::Open(mc_files[iSample].c_str(), "READ");
    // HH: still have the read the individual ROOT file to get the norm histograms
    TH1D* norm = _sampleFile->Get<TH1D>("norm");
    if(!norm){
      MACH3LOG_ERROR("Add a norm KEY to the root file using MakeNormHists.cxx");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    // HH: currently storing the norm_s and POT_s in a map since this is a sample/file specific thing, maybe this could be done in a more elegant way
    norm_map[mc_files[iSample]] = std::vector<double>{norm->GetBinContent(1), pot/norm->GetBinContent(2)};

    // HH: Check whether the file exists, see https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
    int _add_rtn = _data->Add(mc_files[iSample].c_str(), -1);
    if(_add_rtn == 0){
      MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", mc_files[iSample]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    _sampleFile->Close();
  }
  
  
  if(_data){
    // MACH3LOG_INFO("Found \"caf\" tree in {}", mc_files[iSample]);
    MACH3LOG_INFO("Number of entries in TChain: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Failed to create the TChain.");
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

  // now fill the actual variables
  /*
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
  */
  size_t nEntries = static_cast<size_t>(_data->GetEntries());
  dunemcSamples.resize(nEntries);
  _data->GetEntry(0);

  // HH: A map to keep track of negative energies
  std::unordered_map<std::string, int> negative_counts;
  // Initialize the negative counts for each energy variable
  negative_counts["rw_erec_had"] = 0;
  negative_counts["rw_erec_lep"] = 0;
  negative_counts["rw_eRecoN"] = 0;
  negative_counts["rw_eRecoPi0"] = 0;
  negative_counts["rw_sum_ehad"] = 0;
  
  //FILL DUNE STRUCT
  for (unsigned int i = 0; i < nEntries; ++i) { // Loop through tree
    _data->GetEntry(i);

    std::string CurrFileName = _data->GetCurrentFile()->GetName();
    dunemcSamples[i].nupdgUnosc = GetInitPDGFromFileName(CurrFileName);
    dunemcSamples[i].nupdg = GetFinalPDGFromFileName(CurrFileName);
    dunemcSamples[i].OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, dunemcSamples[i].nupdgUnosc, dunemcSamples[i].nupdg));

    // POT stuff
    dunemcSamples[i].norm_s = norm_map[CurrFileName][0]; // Norm in sample
    dunemcSamples[i].pot_s = norm_map[CurrFileName][1]; // POT in sample
    
    dunemcSamples[i].rw_cvnnumu = (_cvnnumu);
    dunemcSamples[i].rw_cvnnue = (_cvnnue);
    dunemcSamples[i].rw_cvnnumu_shifted = (_cvnnumu); 
    dunemcSamples[i].rw_cvnnue_shifted = (_cvnnue);
    if (iselike) {
      dunemcSamples[i].rw_erec = (_erec_nue);
      dunemcSamples[i].rw_erec_shifted = (_erec_nue); 
      dunemcSamples[i].rw_erec_had = (_erec_had_nue);
      dunemcSamples[i].rw_erec_lep = (_erec_lep_nue);
    } else {
      dunemcSamples[i].rw_erec = (_erec); 
      dunemcSamples[i].rw_erec_shifted = (_erec); 
      dunemcSamples[i].rw_erec_had = (_erec_had); 
      dunemcSamples[i].rw_erec_lep = (_erec_lep); 
    }
    
    dunemcSamples[i].rw_eRecoP = (_eRecoP); 
    dunemcSamples[i].rw_eRecoPip = (_eRecoPip); 
    dunemcSamples[i].rw_eRecoPim = (_eRecoPim); 
    dunemcSamples[i].rw_eRecoPi0 = (_eRecoPi0); 
    dunemcSamples[i].rw_eRecoN = (_eRecoN); 
    dunemcSamples[i].rw_LepE = (_LepE); 
    dunemcSamples[i].rw_eP = (_eP); 
    dunemcSamples[i].rw_ePip = (_ePip); 
    dunemcSamples[i].rw_ePim = (_ePim); 
    dunemcSamples[i].rw_ePi0 = (_ePi0); 
    dunemcSamples[i].rw_eN = (_eN);

    // HH: Add checks to make sure the energies are not negative
    if (dunemcSamples[i].rw_erec_had < 0) {
      dunemcSamples[i].rw_erec_had = 0;
      negative_counts["rw_erec_had"]++;
    }
    if (dunemcSamples[i].rw_erec_lep < 0) {
      dunemcSamples[i].rw_erec_lep = 0;
      negative_counts["rw_erec_lep"]++;
    }
    if (dunemcSamples[i].rw_eRecoN < 0) {
      dunemcSamples[i].rw_eRecoN = 0;
      negative_counts["rw_eRecoN"]++;
    }
    if (dunemcSamples[i].rw_eRecoPi0 < 0) {
      dunemcSamples[i].rw_eRecoPi0 = 0;
      negative_counts["rw_eRecoPi0"]++;
    } 

    dunemcSamples[i].rw_erec_had_sqrt = sqrt(dunemcSamples[i].rw_erec_had);
    dunemcSamples[i].rw_erec_lep_sqrt = sqrt(dunemcSamples[i].rw_erec_lep);
    dunemcSamples[i].rw_eRecoN_sqrt = sqrt(dunemcSamples[i].rw_eRecoN);
    dunemcSamples[i].rw_eRecoPi0_sqrt = sqrt(dunemcSamples[i].rw_eRecoPi0);

    dunemcSamples[i].rw_sum_ehad = dunemcSamples[i].rw_eRecoP + dunemcSamples[i].rw_eRecoPip + dunemcSamples[i].rw_eRecoPim;
    if (dunemcSamples[i].rw_sum_ehad < 0) {
      dunemcSamples[i].rw_sum_ehad = 0;
      negative_counts["rw_sum_ehad"]++;
    }
    dunemcSamples[i].rw_sum_ehad_sqrt = sqrt(dunemcSamples[i].rw_sum_ehad);

    dunemcSamples[i].rw_etru = (_ev);
    dunemcSamples[i].rw_isCC = _isCC;
    // dunemcSamples[i].rw_nuPDGunosc = _nuPDGunosc;
    // dunemcSamples[i].rw_nuPDG = _nuPDG;
    dunemcSamples[i].rw_berpaacvwgt = (_BeRPA_cvwgt);
    dunemcSamples[i].rw_vtx_x = (_vtx_x);
    dunemcSamples[i].rw_vtx_y = (_vtx_y);
    dunemcSamples[i].rw_vtx_z = (_vtx_z);

    dunemcSamples[i].rw_trueccnumu = static_cast<double>(dunemcSamples[i].rw_isCC==1 && abs(dunemcSamples[i].nupdg)==14);
    dunemcSamples[i].rw_trueccnue = static_cast<double>(dunemcSamples[i].rw_isCC==1 && abs(dunemcSamples[i].nupdg)==12);
    
    //Assume everything is on Argon40 for now....
    dunemcSamples[i].Target = kTarget_Ar;
    
    int M3Mode = Modes->GetModeFromGenerator(std::abs(_mode));
    if (!_isCC) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    dunemcSamples[i].mode = M3Mode;
    
    dunemcSamples[i].flux_w = 1.0;
  }

  // HH: Give a warning if any negative energies were found
  for (const auto& pair : negative_counts) {
    if (pair.second > 0) {
      MACH3LOG_WARN("Found {} negative values for {} in sample {}", pair.second, pair.first, GetSampleName());
    }
  }
  
  delete _data;
  return static_cast<int>(nEntries);
}

const double* SampleHandlerBeamFD::GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent) {
  double* KinematicValue = nullptr;

  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].rw_etru); 
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].rw_erec_shifted);
    break;
  case kTrueXPos:
    KinematicValue = &(dunemcSamples[iEvent].rw_vtx_x);
    break;
  case kTrueYPos:
    KinematicValue = &(dunemcSamples[iEvent].rw_vtx_y);
    break;
  case kTrueZPos:
    KinematicValue = &(dunemcSamples[iEvent].rw_vtx_z);
    break;
  case kCVNNumu:
    KinematicValue = &(dunemcSamples[iEvent].rw_cvnnumu_shifted);
    break;
  case kCVNNue:
    KinematicValue = &(dunemcSamples[iEvent].rw_cvnnue_shifted);
    break;
  case kM3Mode:
    KinematicValue = &(dunemcSamples[iEvent].mode);
    break;
  case kOscChannel:
    KinematicValue = &(dunemcSamples[iEvent].OscChannelIndex);
    break;
  case kIsFHC:
    KinematicValue = &(isFHC);
    break;
  case kTrueCCnue: 
	KinematicValue = &(dunemcSamples[iEvent].rw_trueccnue);
 	break;
  case kTrueCCnumu: 
	KinematicValue = &(dunemcSamples[iEvent].rw_trueccnumu);
 	break;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* SampleHandlerBeamFD::GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter)); 
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

const double* SampleHandlerBeamFD::GetPointerToKinematicParameter(double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamFD::ReturnKinematicParameter(int KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamFD::ReturnKinematicParameter(std::string KinematicParameter, int iEvent) {
 return *GetPointerToKinematicParameter(KinematicParameter, iEvent);
}

void SampleHandlerBeamFD::SetupFDMC() {
  // dunemc_base *duneobj = &(dunemcSamples[iSample]);
  // FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);  
  
  for(int iEvent = 0 ;iEvent < int(GetNEvents()); ++iEvent) {
    MCSamples[iEvent].rw_etru = &(dunemcSamples[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunemcSamples[iEvent].mode);
    MCSamples[iEvent].Target = &(dunemcSamples[iEvent].Target); 
    MCSamples[iEvent].isNC = !(dunemcSamples[iEvent].rw_isCC);
    MCSamples[iEvent].nupdg = &(dunemcSamples[iEvent].nupdg);
    MCSamples[iEvent].nupdgUnosc = &(dunemcSamples[iEvent].nupdgUnosc);
  }
  
}
 
/*
std::vector<double> SampleHandlerBeamFD::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
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
*/