#include <TROOT.h>

#include "samplePDFDUNEBeamFD.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"

samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_, covarianceXsec* xsec_cov_) : samplePDFFDBase(mc_version_, xsec_cov_) {
  //Call insitialise in samplePDFFD
  Initialise();
}

samplePDFDUNEBeamFD::~samplePDFDUNEBeamFD() {
}

void samplePDFDUNEBeamFD::Init() {
  if (CheckNodeExists(SampleManager->raw(), "DUNESampleBools", "iselike" )) {
	iselike = SampleManager->raw()["DUNESampleBools"]["iselike"].as<bool>();
  } else{
    MACH3LOG_ERROR("Did not find DUNESampleBools:iselike in {}, please add this", SampleManager->GetFileName());
	throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (CheckNodeExists(SampleManager->raw(), "POT")) {
    pot = SampleManager->raw()["POT"].as<double>();
  } else{
    MACH3LOG_ERROR("POT not defined in {}, please add this!", SampleManager->GetFileName());
	throw MaCh3Exception(__FILE__, __LINE__);
  }
 
  tot_escale_fd_pos = -999;
  tot_escale_sqrt_fd_pos = -999;
  tot_escale_invsqrt_fd_pos = -999;
  had_escale_fd_pos = -999;
  had_escale_sqrt_fd_pos = -999;
  had_escale_invsqrt_fd_pos = -999;
  mu_escale_fd_pos = -999;
  mu_escale_sqrt_fd_pos = -999;
  mu_escale_invsqrt_fd_pos = -999;
  n_escale_fd_pos = -999;
  n_escale_sqrt_fd_pos = -999;
  n_escale_invsqrt_fd_pos = -999;
  em_escale_fd_pos = -999;
  em_escale_sqrt_fd_pos = -999;
  em_escale_invsqrt_fd_pos = -999;
  had_res_fd_pos = -999;
  mu_res_fd_pos = -999;
  n_res_fd_pos = -999;
  em_res_fd_pos = -999;
  cvn_numu_fd_pos = -999;
  cvn_nue_fd_pos = -999;

  // create dunemc storage
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }

  nFDDetectorSystPointers = funcParsIndex.size();
  std::unordered_map<std::string, const double*> FDDetectorSystPointersMap;
  FDDetectorSystPointers = std::vector<const double*>(nFDDetectorSystPointers);

  for(auto FuncPar_i  = 0 ; FuncPar_i < funcParsIndex.size() ; ++FuncPar_i){
	FDDetectorSystPointersMap.insert(std::pair<std::string, const double*>(funcParsNames.at(FuncPar_i), XsecCov->retPointer(funcParsIndex.at(FuncPar_i))));
  }

  /*
  int func_it = 0;
  for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++func_it) {
    std::string name = funcParsNames.at(func_it);
    
    if (name == "TotalEScaleFD") {
      tot_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_fd_pos);
    }
    else if (name == "TotalEScaleSqrtFD") {
      tot_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_sqrt_fd_pos);
    }
    else if (name == "TotalEScaleInvSqrtFD") {
      tot_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_invsqrt_fd_pos);
    }
    else if (name == "HadEScaleFD") {
      had_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_fd_pos);
    }
    else if (name == "HadEScaleSqrtFD") {
      had_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_sqrt_fd_pos);
    }
    else if (name == "HadEScaleInvSqrtFD") {
      had_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_invsqrt_fd_pos);
    }
    else if (name == "MuEScaleFD") {
      mu_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_fd_pos);
    }
    else if (name == "MuEScaleSqrtFD") {
      mu_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_sqrt_fd_pos);
    }
    else if (name == "MuEScaleInvSqrtFD") {
      mu_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_invsqrt_fd_pos);
    }
    else if (name == "NEScaleFD") {
      n_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_fd_pos);
    }
    else if (name == "NEScaleSqrtFD") {
      n_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_sqrt_fd_pos);
    }
    else if (name == "NEScaleInvSqrtFD") {
      n_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_invsqrt_fd_pos);
    }
    else if (name == "EMEScaleFD") {
      em_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_fd_pos);
    }
    else if (name == "EMEScaleSqrtFD") {
      em_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_sqrt_fd_pos);
    }
    else if (name == "EMEScaleInvSqrtFD") {
      em_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_invsqrt_fd_pos);
    }
    else if (name == "HadResFD") {
      had_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_res_fd_pos);
    }
    else if (name == "MuResFD") {
      mu_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_res_fd_pos);
    }
    else if (name == "NResFD") {
      n_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_res_fd_pos);
    }
    else if (name == "EMResFD") {
      em_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_res_fd_pos);
    }
    else if (name == "CVNNumuFD") {
      cvn_numu_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(cvn_numu_fd_pos);
    }
    else if (name == "CVNNueFD") {
      cvn_nue_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(cvn_nue_fd_pos);
    }
    
    else { 
      std::cerr << "Found a functional parameter which wasn't specified in the xml | samplePDFDUNEBeamFD:" << name << std::endl;
      throw;
    }
  }
*/
  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void samplePDFDUNEBeamFD::SetupSplines() {

  ///@todo move all of the spline setup into core
  if(XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0){
	MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
	splinesDUNE* DUNESplines = new splinesDUNE(XsecCov);
	splineFile = (splineFDBase*)DUNESplines;
	InitialiseSplineObject();
  }
  else{
	MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
	splineFile = nullptr;
  }

  return;
}

void samplePDFDUNEBeamFD::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j] = new const double*[MCSamples[i].ntotal_weight_pointers[j]];
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
  const char *sampleFile = (mtuple_files[iSample]).c_str();
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  int nutype = sample_nutype[iSample];
  int oscnutype = sample_oscnutype[iSample];
  bool signal = sample_signal[iSample];
  
  std::cout << "-------------------------------------------------------------------" << std::endl;
  MACH3LOG_INFO("input file: {}", sampleFile);
  
  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("caf");

  if(_data){
    MACH3LOG_INFO("Found \"caf\" tree in {}", sampleFile);
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
	MACH3LOG_ERROR("Could not find \"caf\" tree in {}", sampleFile);
	throw MaCh3Exception(__FILE__, __LINE__);
  }
  
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

  TH1D* norm = (TH1D*)_sampleFile->Get("norm");
  if(!norm){
    std::cout<< "Add a norm KEY to the root file using MakeNormHists.cxx"<<std::endl;
    std::cout << "Ignoring for now" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // now fill the actual variables
  duneobj->norm_s = norm->GetBinContent(1);
  duneobj->pot_s = pot/norm->GetBinContent(2);

  std::cout<< "pot_s = " << duneobj->pot_s << std::endl;
  std::cout<< "norm_s = " << duneobj->norm_s << std::endl;

  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;

  std::cout << "signal: " << duneobj->signal << std::endl;
  std::cout << "nevents: " << duneobj->nEvents << std::endl;

  // allocate memory for dunemc variables
  duneobj->rw_cvnnumu.resize(duneobj->nEvents);
  duneobj->rw_cvnnue.resize(duneobj->nEvents);
  duneobj->rw_cvnnumu_shifted.resize(duneobj->nEvents);
  duneobj->rw_cvnnue_shifted.resize(duneobj->nEvents);
  duneobj->rw_etru.resize(duneobj->nEvents);
  duneobj->rw_erec.resize(duneobj->nEvents);
  duneobj->rw_erec_shifted.resize(duneobj->nEvents);
  duneobj->rw_erec_had.resize(duneobj->nEvents);
  duneobj->rw_erec_lep.resize(duneobj->nEvents);

  duneobj->rw_eRecoP.resize(duneobj->nEvents);
  duneobj->rw_eRecoPip.resize(duneobj->nEvents);
  duneobj->rw_eRecoPim.resize(duneobj->nEvents);
  duneobj->rw_eRecoPi0.resize(duneobj->nEvents);
  duneobj->rw_eRecoN.resize(duneobj->nEvents);

  duneobj->rw_LepE.resize(duneobj->nEvents);
  duneobj->rw_eP.resize(duneobj->nEvents);
  duneobj->rw_ePip.resize(duneobj->nEvents);
  duneobj->rw_ePim.resize(duneobj->nEvents);
  duneobj->rw_ePi0.resize(duneobj->nEvents);
  duneobj->rw_eN.resize(duneobj->nEvents);

  duneobj->rw_theta.resize(duneobj->nEvents);
  duneobj->flux_w.resize(duneobj->nEvents);
  duneobj->rw_isCC.resize(duneobj->nEvents);
  duneobj->rw_nuPDGunosc.resize(duneobj->nEvents);
  duneobj->rw_nuPDG.resize(duneobj->nEvents);
  duneobj->rw_berpaacvwgt.resize(duneobj->nEvents); 
  duneobj->rw_vtx_x.resize(duneobj->nEvents);
  duneobj->rw_vtx_y.resize(duneobj->nEvents);
  duneobj->rw_vtx_z.resize(duneobj->nEvents);

  duneobj->global_bin_number.resize(duneobj->nEvents);

  duneobj->mode.resize(duneobj->nEvents);
  duneobj->Target.resize(duneobj->nEvents);

  _data->GetEntry(0);

  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    _data->GetEntry(i);
    duneobj->rw_cvnnumu[i] = _cvnnumu;
    duneobj->rw_cvnnue[i] = _cvnnue;
    duneobj->rw_cvnnumu_shifted[i] = _cvnnumu; 
    duneobj->rw_cvnnue_shifted[i] = _cvnnue;
    if (iselike) {
      duneobj->rw_erec[i] = _erec_nue;
      duneobj->rw_erec_shifted[i] = _erec_nue; 
      duneobj->rw_erec_had[i] = _erec_had_nue;
      duneobj->rw_erec_lep[i] = _erec_lep_nue;
    } else {
      duneobj->rw_erec[i] = _erec; 
      duneobj->rw_erec_shifted[i] = _erec; 
      duneobj->rw_erec_had[i] = _erec_had; 
      duneobj->rw_erec_lep[i] = _erec_lep; 
    }
    
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
    
    duneobj->rw_etru[i] = _ev;
    duneobj->rw_theta[i] = _LepNuAngle;
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;
    duneobj->rw_vtx_x[i] = _vtx_x;
    duneobj->rw_vtx_y[i] = _vtx_y;
    duneobj->rw_vtx_z[i] = _vtx_z;

    duneobj->global_bin_number[i] = GetGenericBinningGlobalBinNumber(iSample, i);
    
    //Assume everything is on Argon for now....
    duneobj->Target[i] = 40;
    
    int mode= TMath::Abs(_mode);       
    duneobj->mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
    
    duneobj->flux_w[i] = 1.0;
  }
  
  _sampleFile->Close();
  return duneobj->nEvents;
}

TH1D* samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill!=-1) {
    if (kChannelToFill>dunemcSamples.size()) {
      std::cout << "Required channel is not available. kChannelToFill should be between 0 and " << dunemcSamples.size() << std::endl;
      std::cout << "kChannelToFill given:" << kChannelToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill!=-1) {
    if (kModeToFill>kMaCh3_nModes) {
      std::cout << "Required mode is not available. kModeToFill should be between 0 and " << kMaCh3_nModes << std::endl;
      std::cout << "kModeToFill given:" << kModeToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw;
    }
    fMode = true;
  } else {
    fMode = false;
  }

  std::vector< std::vector<double> > SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = kM3Mode;
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill+1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = kOscChannel;
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill+1;
    SelectionVec.push_back(SelecChannel);
  }

  return get1DVarHist(Var1,SelectionVec,WeightStyle,Axis);
}

/*! DB New version of get1DVarHist which only fills histogram with events passing IsEventSelected
 * This works by having the Selection vector, where each component of Selection is a 2 or 3 length vector
 * If Selection[i].size()==3, Selection[i][0] is the ND280KinematicType which is being cut, and only events with ND280KinematicType values between Selection[i][1] and Selection[i][2] are accepted
 */
TH1D* samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1,std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis) {

  Selection = SelectionVec;

  for (unsigned int iStoredSelection=0;iStoredSelection<StoredSelection.size();iStoredSelection++) {
    Selection.push_back(StoredSelection[iStoredSelection]);
  }

  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
    if (Selection[iSelection].size()!=3) {
      std::cerr << "Selection Vector[" << iSelection << "] is not formed correctly. Expect size == 3, given:" << Selection[iSelection].size() << std::endl;
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  //DB Cut on OscChannel in this function due to speed increase from considering duneSamples structure (ie. Array of length NChannels)
  bool fChannel = false;
  int kChannelToFill = -1;
  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
    if (Selection[iSelection][0] == kOscChannel) {
      fChannel = true;
      kChannelToFill = Selection[iSelection][1];
    }
  }

  if (fChannel && kChannelToFill>dunemcSamples.size()) {
    std::cout << "Required channel is not available. kChannelToFill should be between 0 and " << dunemcSamples.size() << std::endl;
    std::cout << "kChannelToFill given:" << kChannelToFill << std::endl;
    std::cout << "Exitting.." << std::endl;
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TH1D* _h1DVar;
  std::vector<double> xBinEdges = ReturnKinematicParameterBinning(Var1);
  _h1DVar = new TH1D("", "", xBinEdges.size()-1, xBinEdges.data());

  //This should be the same as FillArray in core basically, except that
  //events will end up in different bins
  for (int i=0;i<dunemcSamples.size();i++) {
    if (fChannel && (i!=kChannelToFill)) {
      continue;
    }
    for(int j=0;j<dunemcSamples[i].nEvents;j++) {

      //DB Determine which events pass selection
      if (!IsEventSelected(i,j)) {
		continue;
      }

      double Weight = GetEventWeight(i,j);
	  if (WeightStyle==1) {
	    Weight = *(MCSamples[i].osc_w_pointer[j]) * dunemcSamples[i].pot_s * dunemcSamples[i].norm_s * dunemcSamples[i].flux_w[j];
	  }

	  //ETA - not sure about this
	  if (MCSamples[i].xsec_w[j] == 0.) continue;

	  double Var1_Val;

	  Var1_Val = ReturnKinematicParameter(Var1,i,j);
	  if (Var1_Val!=_DEFAULT_RETURN_VAL_) {
		_h1DVar->Fill(Var1_Val,Weight);
	  }
    }
  }

  /* DB: This is commented out be default
  // This code shifts the histogram meaning to Events/Bin Width but this affects the overall integral of the histogram so it should not be used anywhere we care about event rates
  // We could use Hist->Integral("width") but it would require a lot of modification throughout the code

  if (Var1!=kPDFBinning) {
    //_h1DVar->SetBinContent(1,_h1DVar->GetBinContent(0)+_h1DVar->GetBinContent(1));
    //_h1DVar->SetBinContent(_h1DVar->GetNbinsX(),_h1DVar->GetBinContent(_h1DVar->GetNbinsX())+_h1DVar->GetBinContent(_h1DVar->GetNbinsX()+1));

    for (int x=1;x<=_h1DVar->GetNbinsX();x++) {
      _h1DVar->SetBinContent(x,_h1DVar->GetBinContent(x)/_h1DVar->GetXaxis()->GetBinWidth(x));
    }

    _h1DVar->GetYaxis()->SetTitle("Events/Bin Width");
  }
  */

  return _h1DVar;
}

double const& samplePDFDUNEBeamFD::ReturnKinematicParameterByReference(int KinematicParameter, int iSample, int iEvent) {

  switch(KinematicParameter){
  case kTrueNeutrinoEnergy:
    return dunemcSamples[iSample].rw_etru[iEvent]; 
  case kRecoNeutrinoEnergy:
    return dunemcSamples[iSample].rw_erec_shifted[iEvent];
  case kTrueXPos:
    return dunemcSamples[iSample].rw_vtx_x[iEvent];
  case kTrueYPos:
    return dunemcSamples[iSample].rw_vtx_y[iEvent];
  case kTrueZPos:
    return dunemcSamples[iSample].rw_vtx_z[iEvent];
  case kCVNNumu:
    return dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
  case kCVNNue:
    return dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
  case kGlobalBinNumber:
    return dunemcSamples[iSample].global_bin_number[iEvent];
  default:
    std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
    throw;
  }
  }

double samplePDFDUNEBeamFD::ReturnKinematicParameter(int KinematicParameter, int iSample, int iEvent) {
  return ReturnKinematicParameterByReference(KinematicParameter, iSample, iEvent);
}

int samplePDFDUNEBeamFD::ReturnKinematicParameterFromString(std::string KinematicParameterStr){

  if (KinematicParameterStr.find("TrueNeutrinoEnergy") != std::string::npos) {return kTrueNeutrinoEnergy;}
  if (KinematicParameterStr.find("RecoNeutrinoEnergy") != std::string::npos) {return kRecoNeutrinoEnergy;}
  if (KinematicParameterStr.find("TrueXPos") != std::string::npos) {return kTrueXPos;}
  if (KinematicParameterStr.find("TrueYPos") != std::string::npos) {return kTrueYPos;}
  if (KinematicParameterStr.find("TrueZPos") != std::string::npos) {return kTrueZPos;}
  if (KinematicParameterStr.find("CVNNumu") != std::string::npos) {return kCVNNumu;}
  if (KinematicParameterStr.find("CVNNue") != std::string::npos) {return kCVNNue;}
  if (KinematicParameterStr.find("M3Mode") != std::string::npos) {return kM3Mode;}
  if (KinematicParameterStr.find("global_bin_number") != std::string::npos) {return kGlobalBinNumber;}


}

inline std::string samplePDFDUNEBeamFD::ReturnStringFromKinematicParameter(int KinematicParameter) {
  std::string KinematicString = "";
 
  switch(KinematicParameter){
   case kRecoNeutrinoEnergy:
     KinematicString = "RecoNeutrinoEnergy";
	 break;
   case kTrueNeutrinoEnergy:
     KinematicString = "RecoNeutrinoEnergy";
	 break;
   case kTrueXPos:
	 KinematicString= "TrueXPos";
	 break;
   case kTrueYPos:
	 KinematicString= "TrueYPos";
	 break;
   case kTrueZPos:
	 KinematicString= "TrueZPos";
	 break;
   case kCVNNumu:
	 KinematicString = "CVNNumu";
	 break;
   case kCVNNue:
	 KinematicString = "CVNNue";
	 break;
   case kM3Mode:
   KinematicString = "M3Mode";
   break;
   case kGlobalBinNumber:
   KinematicString = "global_bin_number";
   break;
   default:
    break;
  }

  return KinematicString;
}

void samplePDFDUNEBeamFD::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  fdmc_base *fdobj = &(MCSamples[iSample]);  
  
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
  }
}
 
void samplePDFDUNEBeamFD::applyShifts(int iSample, int iEvent) {
   
  //ETA - this is pretty horrific... we need to think of a nicer way to do this.
  //Don't want to add in hard checks on which systematics are defined but also don't want to hard-code
  //the order in which the systematics are specified. All of these functions should have access to the 
  //dunemc struct so they only need to have iSample and iEvent passed to them. Can probably loop over
  //a vector of std::function objects and pass each of them iSample and iEvent.
  /*
   // reset erec back to original value
  dunemcSamples[iSample].rw_erec_shifted[iEvent] = dunemcSamples[iSample].rw_erec[iEvent];

  // reset cvnnumu back to original value
  dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnumu[iEvent];

  // reset cvnnue back to original value
  dunemcSamples[iSample].rw_cvnnue_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnue[iEvent];

  //Calculate values needed
  double sqrtErecHad =  sqrt(dunemcSamples[iSample].rw_erec_had[iEvent]);
  double sqrtErecLep =  sqrt(dunemcSamples[iSample].rw_erec_lep[iEvent]);
  double sqrteRecoPi0 = sqrt(dunemcSamples[iSample].rw_eRecoPi0[iEvent]);
  double sqrteRecoN = sqrt(dunemcSamples[iSample].rw_eRecoN[iEvent]);
  double sumEhad = dunemcSamples[iSample].rw_eRecoP[iEvent] + dunemcSamples[iSample].rw_eRecoPip[iEvent] + dunemcSamples[iSample].rw_eRecoPim[iEvent];
  double sqrtSumEhad = sqrt(sumEhad);

  double invSqrtErecHad =  1/(sqrtErecHad+0.1);
  double invSqrtErecLep =  1/(sqrtErecLep+0.1);
  double invSqrteRecoPi0 =  1/(sqrteRecoPi0+0.1);
  double invSqrteRecoN =  1/(sqrteRecoN+0.1);
  double invSqrtSumEhad =  1/(sqrtSumEhad+0.1);

  bool CCnumu {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14) && dunemcSamples[iSample].nutype==2};
  bool CCnue {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==12) && dunemcSamples[iSample].nutype==1};
  bool NotCCnumu {!(dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14)) && dunemcSamples[iSample].nutype==2};


  TotalEScaleFD(FDDetectorSystPointers[0], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], NotCCnumu);

  TotalEScaleSqrtFD(FDDetectorSystPointers[1], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecHad, sqrtErecLep, NotCCnumu);

  TotalEScaleInvSqrtFD(FDDetectorSystPointers[2], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecHad, invSqrtErecLep, NotCCnumu);

  HadEScaleFD(FDDetectorSystPointers[3], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad);

  HadEScaleSqrtFD(FDDetectorSystPointers[4], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, sqrtSumEhad);

  HadEScaleInvSqrtFD(FDDetectorSystPointers[5], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, invSqrtSumEhad);

  MuEScaleFD(FDDetectorSystPointers[6], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], CCnumu);

  MuEScaleSqrtFD(FDDetectorSystPointers[7], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, CCnumu);

  MuEScaleInvSqrtFD(FDDetectorSystPointers[8], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, CCnumu);

  NEScaleFD(FDDetectorSystPointers[9], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent]);

  NEScaleSqrtFD(FDDetectorSystPointers[10], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], sqrteRecoN);

  NEScaleInvSqrtFD(FDDetectorSystPointers[11], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], invSqrteRecoN);

  EMEScaleFD(FDDetectorSystPointers[12], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], CCnue);

  EMEScaleSqrtFD(FDDetectorSystPointers[13], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, sqrteRecoPi0, CCnue);

  EMEScaleInvSqrtFD(FDDetectorSystPointers[14], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, invSqrteRecoPi0, CCnue);

  HadResFD(FDDetectorSystPointers[15], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoP[iEvent], dunemcSamples[iSample].rw_eRecoPip[iEvent], dunemcSamples[iSample].rw_eRecoPim[iEvent], dunemcSamples[iSample].rw_eP[iEvent], dunemcSamples[iSample].rw_ePip[iEvent], dunemcSamples[iSample].rw_ePim[iEvent]);

  MuResFD(FDDetectorSystPointers[16], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], dunemcSamples[iSample].rw_LepE[iEvent], CCnumu);

  NResFD(FDDetectorSystPointers[17], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], dunemcSamples[iSample].rw_eN[iEvent]);

  EMResFD(FDDetectorSystPointers[18], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_ePi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], dunemcSamples[iSample].rw_LepE[iEvent], CCnue);

  CVNNumuFD(FDDetectorSystPointers[19], &dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent]);

  CVNNueFD(FDDetectorSystPointers[20], &dunemcSamples[iSample].rw_cvnnue_shifted[iEvent]);
  */
}

std::vector<double> samplePDFDUNEBeamFD::ReturnKinematicParameterBinning(int KinematicParameter) {
  std::vector<double> binningVector;

  int nBins = 0;
  double bin_width = 0;
  switch(KinematicParameter){
	case(kRecoNeutrinoEnergy):
	  nBins = 20; 
	  bin_width = 0.5; //GeV
	  break;
	case(kTrueNeutrinoEnergy):
	  nBins = 20; 
	  bin_width = 0.5; //GeV
	  break;
	default:
	  nBins = 10;
	  bin_width = 1.0;
	  break;
  }

  for(int bin_i = 0 ; bin_i < nBins ; bin_i++){
	binningVector.push_back(bin_i*bin_width);
  }

  return binningVector;
}
