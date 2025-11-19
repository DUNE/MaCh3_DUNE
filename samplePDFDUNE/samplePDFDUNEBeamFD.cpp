#include "samplePDFDUNEBeamFD.h"
#include "TVector3.h"
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <TH2.h>

#include <yaml-cpp/yaml.h>
#include <set>
#include <iostream>

#include <fstream>



#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic pop

samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  OscCov = nullptr;
  
  
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
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
}

// === HH: Functional parameters ===
// void samplePDFDUNEBeamFD::TotalEScaleND(const double * par, std::size_t iSample, std::size_t iEvent) {
//   // Total energy scale uncertainties for anything but CC Numu, see:
//   // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
//   //MACH3LOG_INFO("TotalEScaleND par = %g, shift = %g", *par, (*par) * dunemcSamples[iSample].rw_erec_had[iEvent]);
//   std::stringstream ss;
//   ss << "TotalEScaleND par = " << *par
//    << ", shift = " << ((*par) * dunemcSamples[iSample].rw_erec_had[iEvent]);

//   MACH3LOG_INFO(ss.str());

//   std::cout<< "TotalEScaleND =  "  << " dunemcSamples[iSample].rw_erec_had[iEvent] =. " <<  dunemcSamples[iSample].rw_erec_had[iEvent] <<std::endl;
//   dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_had[iEvent];
// }


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
double samplePDFDUNEBeamFD::CalculatePOT() {
  TChain calc_pot_chain("meta");  // Use correct tree name

  std::string pot_branch = "pot";  // Use the correct branch name

  for (int i = 0; i < static_cast<int> (dunemcSamples.size()); ++i) {
      calc_pot_chain.AddFile((mc_files[i]).c_str());
  }

  // Check if the branch exists before proceeding
  if (!calc_pot_chain.GetBranch(pot_branch.c_str())) {
      std::cerr << "Error: Branch " << pot_branch << " not found in the tree!" << std::endl;
      return 0.0;
  }

  double pot_value = 0.0;
  calc_pot_chain.SetBranchAddress(pot_branch.c_str(), &pot_value);

  double sum_pot = 0.0;
  Long64_t nEntries = calc_pot_chain.GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
      calc_pot_chain.GetEntry(i);
      sum_pot += pot_value;
  }

  std::cout << "Summed POT: " << sum_pot << std::endl;
  return sum_pot;
}

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
  TTree* _data = _sampleFile->Get<TTree>("cafTree");
  TTree* _meta = _sampleFile->Get<TTree>("meta");
  if(_data){
    MACH3LOG_INFO("Found \"cafTree\" tree in {}", mc_files[iSample]);
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find \"caf\" tree in {}", mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  _meta = static_cast<TTree*>( _sampleFile->Get("meta"));
  if (!_meta ){
      MACH3LOG_ERROR("Branch 'meta' not found in file: {}", mc_files[iSample]);
      throw MaCh3Exception(__FILE__, __LINE__);
  } 
 
  double newpot = CalculatePOT();
  std::cout<< "CalculatePOT() = " << newpot << " -----------------------------------------------------------------------------" << std::endl;

   
  //Reco Variables
  double _erec;
  double _erec_nue;
  double _erec_had;
  double _erec_had_nue;
  //double _erec_lep;
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
  bool _isNC;
  int _nuPDGunosc;
  int _nuPDG;
  double _LepMomX;
  double _LepMomY;
  double _LepMomZ;
  double _NuMomX;
  double _NuMomY;
  double _NuMomZ;

  double _LepNuAngle;
  double _Elep_reco;

  double _eOther;
  Int_t _nipi0;
  
  double _production_pot = 0.0;  // Explicitly initialize
  double gen_pot = 0.0; //set the sum of the pot from each file to be 0,befor any are read in

  if(_meta){
    std::cout<<"Found meta tree" << std::endl;

    _meta->SetBranchStatus("*", 0);
    _meta->SetBranchStatus("pot", 1);
    _meta->SetBranchAddress("pot", &_production_pot);

    for (int i = 0; i < _meta->GetEntries(); i++) {
      _meta->GetEntry(i);
      if (_production_pot > 1e30) {  // Arbitrary threshold for large numbers
        std::cerr << "WARNING: Unusually large pot value at entry " << i << ": " 
                  << _production_pot << std::endl;
      }
      gen_pot = gen_pot + _production_pot;
  }
  }
  else{
    MACH3LOG_ERROR("Could not find \"meta\" tree in {}", mc_files[iSample]);
    gen_pot =  3.85e21;
    std::cout << "instead setting gen_pot to be = " << gen_pot << std::endl;
  }
  std::cout << " final gen_pot  = " << gen_pot << std::endl;


  Int_t _reco_numu;
  Int_t _muon_contained;
  Int_t _muon_tracker;
  Double_t _Ehad_veto;


  
  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
 // _data->SetBranchStatus("Ev_reco_numu", 1);
 // _data->SetBranchAddress("Ev_reco_numu", &_erec);
  _data->SetBranchStatus("Ev_reco_nue", 1);
  _data->SetBranchAddress("Ev_reco_nue", &_erec_nue);
  _data->SetBranchStatus("Ev_reco", 1);
  _data->SetBranchAddress("Ev_reco", &_erec);
  _data->SetBranchStatus("Elep_reco", 1);
  _data->SetBranchAddress("Elep_reco", &_Elep_reco);


  _data->SetBranchStatus("RecoHadEnNumu", 1);
  _data->SetBranchAddress("RecoHadEnNumu", &_erec_had);
  _data->SetBranchStatus("RecoHadEnNue", 1);
  _data->SetBranchAddress("RecoHadEnNue", &_erec_had_nue);
  //_data->SetBranchStatus("RecoLepEnNumu", 1);
  //_data->SetBranchAddress("RecoLepEnNumu", &_erec_lep);
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

 ///////////////////////////

  _data->SetBranchStatus("reco_numu", 1);
  _data->SetBranchAddress("reco_numu", &_reco_numu);

  _data->SetBranchStatus("muon_contained", 1);
  _data->SetBranchAddress("muon_contained", &_muon_contained);

  _data->SetBranchStatus("muon_tracker", 1);
  _data->SetBranchAddress("muon_tracker", &_muon_tracker);

  _data->SetBranchStatus("Ehad_veto", 1);
  _data->SetBranchAddress("Ehad_veto", &_Ehad_veto);

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

  _data->SetBranchStatus("eOther", 1);
  _data->SetBranchAddress("eOther", &_eOther);   //nipi0 

  _data->SetBranchStatus("nipi0", 1);
  _data->SetBranchAddress("nipi0", &_nipi0);

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


  _data->SetBranchStatus("LepMomX", 1);
  _data->SetBranchAddress("LepMomX", &_LepMomX);  
  _data->SetBranchStatus("LepMomY", 1);
  _data->SetBranchAddress("LepMomY", &_LepMomY);  
  _data->SetBranchStatus("LepMomZ", 1);
  _data->SetBranchAddress("LepMomZ", &_LepMomZ);  
  _data->SetBranchStatus("LepMomZ", 1);
  _data->SetBranchAddress("LepMomZ", &_LepMomZ);  

  _data->SetBranchStatus("NuMomX", 1);
  _data->SetBranchAddress("NuMomX", &_NuMomX);  
  _data->SetBranchStatus("NuMomY", 1);
  _data->SetBranchAddress("NuMomY", &_NuMomY);  
  _data->SetBranchStatus("NuMomZ", 1);
  _data->SetBranchAddress("NuMomZ", &_NuMomZ);  
  _data->SetBranchStatus("NuMomZ", 1);
  _data->SetBranchAddress("NuMomZ", &_NuMomZ);  
  _data->SetBranchStatus("LepNuAngle", 1);
  _data->SetBranchAddress("LepNuAngle", &_LepNuAngle);

  _data->SetBranchStatus("isNC", 1);
  _data->SetBranchAddress("isNC", &_isNC);

  

  TH1D* norm = _sampleFile->Get<TH1D>("norm");
  if(!norm){
    MACH3LOG_ERROR("Add a norm KEY to the root file using MakeNormHists.cxx");
    //throw MaCh3Exception(__FILE__, __LINE__);
    norm = new TH1D("norm","",1,0,1);
    norm->SetBinContent(1,1);
    duneobj->norm_s = 1.0; 
    duneobj->pot_s = (pot)/(newpot);
    std::cout << "(pot)/(newpot)" << (pot)/(newpot) << std::endl;
  }
  else{
    duneobj->norm_s = norm->GetBinContent(1);
    duneobj->pot_s = pot/norm->GetBinContent(2);
  }
  std::cout << "pot = " << (pot)<< std::endl;
  std::cout << "pot_s = " << duneobj->pot_s << std::endl;

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

  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new double[duneobj->nEvents];
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

  duneobj->true_q0 = new double[duneobj->nEvents];
  duneobj->true_q3 = new double[duneobj->nEvents];
  duneobj->rw_pt = new double[duneobj->nEvents];
  duneobj->rw_pz = new double[duneobj->nEvents];

  duneobj->global_bin_number = new double[duneobj->nEvents];
  duneobj->template_global_bin_number = new double[duneobj->nEvents];
  duneobj->p_lep = new double[duneobj->nEvents];
  duneobj->theta_lep = new double[duneobj->nEvents];
  duneobj->rw_yrec = new double[duneobj->nEvents];
  
  duneobj->enu_proxy_minus_enutrue = new double[duneobj->nEvents];
  duneobj->ERec_QE = new double[duneobj->nEvents];
  duneobj->erec_proxy = new double[duneobj->nEvents];
  duneobj->erec_proxy_minus_enu = new double[duneobj->nEvents];
  duneobj->eHad_av = new double[duneobj->nEvents];///////////////////////new variable :) bit like the truth one in GENIE

  duneobj->reco_numu = new bool[duneobj->nEvents];
  duneobj->muon_contained = new bool[duneobj->nEvents];
  duneobj->muon_tracker = new bool[duneobj->nEvents];
  duneobj->Ehad_veto = new float[duneobj->nEvents];
  duneobj->isNC = new bool[duneobj->nEvents];
  duneobj->enurec_minus_enutrue = new double[duneobj->nEvents];
  duneobj->relative_enu_bias = new double[duneobj->nEvents];
  duneobj->enu_bias = new double[duneobj->nEvents];

  _data->GetEntry(0);
  bool need_global_bin_numbers = (XVarStr == "global_bin_number");


  TH2D* globalMap = new TH2D("global_bin_map",
                           "Global vs Template Global;Global bin;Template bin",
                           576, 0.5, 576.5,
                           480, 0.5, 480.5);

  // Open file at the start of your program
  std::ofstream outfile("bin_numbers.txt");
  if(!outfile.is_open()){
      std::cerr << "Cannot open output file!" << std::endl;
  }


  //FILL DUNE STRUCT

  int conditionCounter = 0;
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    
    _data->GetEntry(i);

    duneobj->reco_numu[i] = _reco_numu;
    duneobj->muon_contained[i] = _muon_contained;
    duneobj->muon_tracker[i] = _muon_tracker;
    duneobj->Ehad_veto[i] = static_cast<float>(_Ehad_veto);

    duneobj->nupdg[i] = sample_nupdg[iSample];
    duneobj->nupdgUnosc[i] = sample_nupdgunosc[iSample];    
    
    duneobj->rw_cvnnumu[i] = (_cvnnumu);
    duneobj->rw_cvnnue[i] = (_cvnnue);
    duneobj->rw_cvnnumu_shifted[i] = (_cvnnumu); 
    duneobj->rw_cvnnue_shifted[i] = (_cvnnue);
    double el = 0.0;

    duneobj->rw_erec[i] = (_erec); //------------------------------------original version!!
    
    duneobj->rw_erec_shifted[i] = (_erec); 
    //std::cout << "rw_erec_shifted[i]  = " << (_erec) << std::endl;
    //std::cout << "rw_erec_had[i] = " <<   (_erec - _Elep_reco) << std::endl;
    
    duneobj->rw_yrec[i] = ((_erec -_Elep_reco)/_erec);

    //std::cout<< " rw_yrec"  << ((_erec -_Elep_reco)/_erec) << std::endl;
    //std::cout<< " erec  = " <<_erec << std::endl;
    //std::cout<< " elep_reco  = " << _Elep_reco << std::endl;
    duneobj->enurec_minus_enutrue[i] = (_erec) - (_ev);
    //if (iselike) {
      //duneobj->rw_erec[i] = (_erec);
      //duneobj->rw_erec_shifted[i] = (_erec_nue); 
      //duneobj->rw_erec_had[i] = (_erec -_Elep_reco);
      
      //duneobj->rw_erec_lep[i] = (_Elep_reco);
      //duneobj->enurec_minus_enutrue[i] = (_erec) - (_ev);
      //el = (_erec_lep_nue);

      //duneobj->erec_proxy[i] = (_erec_lep_nue) + (_erec_had_nue);
      //duneobj->erec_proxy_minus_enu[i] = (_erec_lep_nue) + (_erec_had_nue) - (_ev);
      //std::cout << "erec_proxy_minus_enu = " <<   _Elep_reco  - (_ev) << std::endl;
      //std::cout << "erec_proxy_minus_enu = " <<  (_erec_lep) + (_erec_had) - (_ev) << std::endl;
      //std::cout << "erec_lep = " << _Elep_reco << std::endl;
      //std::cout << "erec_had = " <<   (_erec -_Elep_reco) << std::endl;

    //} else {
      //duneobj->rw_erec[i] = (_erec); 
      //duneobj->rw_erec_shifted[i] = (_erec); 
      //duneobj->rw_erec_had[i] = (_erec -_Elep_reco); 
      //duneobj->rw_erec_lep[i] = (_Elep_reco); 
      
      //duneobj->enurec_minus_enutrue[i] = (_erec) - (_ev);
      //el = (_erec_lep);

      //duneobj->erec_proxy[i] = (_erec_lep) + (_erec_had);
    //}

    //////////////////////////////////Also EHadav which is a replacment for not being able to access the GENIE tree in the OA CAF's 
    ///////////EHadAv = eP + ePip + ePim +ePi0 + eOther + (nipi0*0.1349)
    double eHad_truth =  _eP + _ePip + _ePim + _ePi0 + _eOther + (_nipi0 *0.1349); 
    duneobj->eHad_av[i]= eHad_truth;
    //std::cout << "eHad_truth = " << eHad_truth << std::endl;
    //std::cout << "eHad_av = " << _eP + _ePip + _ePim + _ePi0 + _eOther + (_nipi0 *0.1349) << std::endl;
    //duneobj->rw_erec[i] = _Elep_reco + eHad_truth;//////////////////////new test for stephen and laura :)

    

    //std::cout << "Original Enu rec = rw_erec = " <<  (_erec) << std::endl;
    //std::cout << "Enu rec = _Elep_reco + _erec_had = " << _Elep_reco << " + " << eHad_truth << " = " << (_Elep_reco + eHad_truth) << std::endl;

    


   
    duneobj->rw_erec_had[i] = (_erec - _Elep_reco);
    duneobj->enu_proxy_minus_enutrue[i] = (_LepE) + eHad_truth - _ev;
    //std::cout << "         _LepE  =  "  << _LepE << std::endl;
    //std::cout<< "_Elep_rec = " << (_Elep_reco) << std::endl;
    //std::cout<< "rw_erec_had[i] = " <<  (_erec -_Elep_reco) << std::endl;
    //std::cout<< "enu true[i] = " << ( _ev) << std::endl;
    //std::cout<< "yrec[i] = " << ((_erec -_Elep_reco)/_erec) << std::endl;
    //std::cout<< "erec_proxy_minus_enu[i] = " <<  (_Elep_reco) + eHad_truth - _ev << std::endl;
    duneobj->isNC[i] = _isNC;
   
    /*
    if(_Elep_reco > 1.0 * _erec ){
      //conditionCounter++;
      std::cout << "condition _Elep_reco > 1.0 * _erec is satisfied..." << std::endl; 
      std::cout<< "Original Enu. rec = " << _erec << std::endl;
      std::cout<< "Erec lep.  = " << _Elep_reco << std::endl;
      std::cout<< "Ehad av (truth Ehad = )  = " << eHad_truth << std::endl;
    }
    */
    if(_Elep_reco != 0.0 ){
       duneobj->rw_erec_lep[i] = (_Elep_reco);
    }
    
    
    duneobj->rw_eRecoP[i] = (_eRecoP); 
    duneobj->rw_eRecoPip[i] = (_eRecoPip); 
    duneobj->rw_eRecoPim[i] = (_eRecoPim); 
    duneobj->rw_eRecoPi0[i] = (_eRecoPi0); 
    duneobj->rw_eRecoN[i] = (_eRecoN); 
    
    duneobj->rw_LepE[i] =(_LepE); 
    duneobj->rw_eP[i] = (_eP); 
    duneobj->rw_ePip[i] = (_ePip); 
    duneobj->rw_ePim[i] = (_ePim); 
    duneobj->rw_ePi0[i] = (_ePi0); 
    duneobj->rw_eN[i] = (_eN); 
    
    duneobj->true_q0[i] = (_ev - _LepE);
    TVector3 nuMom(_NuMomX, _NuMomY, _NuMomZ);
    TVector3 nuMomNorm = nuMom.Unit(); // Normalized vector
    duneobj->true_q3[i] = (TVector3{_NuMomX, _NuMomY, _NuMomZ} -TVector3{_LepMomX, _LepMomY, _LepMomZ}).Mag();
    duneobj->rw_etru[i] = (_ev);
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = (_BeRPA_cvwgt);
    duneobj->rw_vtx_x[i] = (_vtx_x);
    duneobj->rw_vtx_y[i] = (_vtx_y);
    duneobj->rw_vtx_z[i] = (_vtx_z);

    duneobj->rw_pt[i] = TVector3(_LepMomX, _LepMomY, _LepMomZ).Dot(nuMomNorm);
    duneobj->rw_pz[i] =  (TVector3(_LepMomX, _LepMomY, _LepMomZ).Cross(nuMomNorm)).Mag();
    
    //Assume everything is on Argon40 for now....
    duneobj->Target[i] = kTarget_Ar;

    duneobj->theta_lep[i] = (_LepNuAngle);
    duneobj->p_lep[i] =(TVector3{_LepMomX, _LepMomY, _LepMomZ}).Mag();

    duneobj->relative_enu_bias[i] = ((_LepE) + eHad_truth - _ev)/(_ev);
    duneobj->enu_bias[i] = ((_LepE) + eHad_truth - _ev);
    

    //Longer calculation for ERecQE-------------------------------------------------------------------------------
    constexpr double V = 0;        // 0 binding energy for now
    constexpr double mn = 939.565; // neutron mass
    constexpr double mp = 938.272; // proton mass
    double mN_eff = mn - V;
    double mN_oth = mp;

      if (_nuPDGunosc < 0) { // if anti-neutrino, swap target/out masses
        mN_eff = mp - V;
        mN_oth = mn;
      }
     
      // this is funky, but don't be scared, it defines an annonymous function
      // in place that grabs the lepton mass in MeV when given the neutrino PDG
      // and whether the interaction was CC or NC and then immediately calls it.
      // It's basically a generalisation of the ternary operator.
      double ml =
        [](int nupdg, bool isCC) {
          switch (std::abs(nupdg)) {
          case 12: return isCC ? 0.511 : 0;
          case 14: return isCC ? 105.66 : 0;
          case 16: return isCC ? 1777.0 : 0;
          default:
            std::cerr << "Warning: Unexpected PDG code " << nupdg << " passed to ml lambda.\n";
            assert(false && "Unexpected neutrino PDG code in ml lambda");
            return 0.0;
          }
        }(_nuPDGunosc, _isCC);

      double pl = std::sqrt(el*el - ml*ml); // momentum of lepton

      double rEnu =
          (2 * mN_eff * el - ml * ml + mN_oth * mN_oth - mN_eff * mN_eff) /
          (2 * (mN_eff - el +
                pl * std::cos((_LepNuAngle))));

    duneobj->ERec_QE[i] = rEnu;
    
    //------------------------------------------------------------------------------------------------------------
    int M3Mode = Modes->GetModeFromGenerator(std::abs(_mode));
    if (!_isCC) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    duneobj->mode[i] = M3Mode;
    
    duneobj->flux_w[i] = 1.0;

      if(need_global_bin_numbers){
        //Observable binning
        std::vector<double> global_enu_bins = {0.0000, 0.1364, 0.2727, 0.4091, 0.5455, 0.6818, 0.8182, 0.9545, 1.0909, 1.2273, 1.3636, 1.5000, 1.6364, 1.7727, 1.9091, 2.0455, 2.1818, 2.3182,2.4545, 2.5909, 2.7273, 2.8636, 3.0000, 5, 6};
        std::vector<double> global_elep_bins = {0.0000, 0.1364, 0.2727, 0.4091, 0.5455, 0.6818, 0.8182, 0.9545, 1.0909, 1.2273, 1.3636, 1.5000, 1.6364, 1.7727, 1.9091, 2.0455, 2.1818, 2.3182,2.4545, 2.5909, 2.7273, 2.8636, 3.0000, 5, 6};
        
        //Template parameter binning
        std::vector<double> tmplt_enu_bins = {0.0, 1.0, 1.25, 1.5,1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 5.0, 6.0 ,10.0};
        std::vector<double> tmplt_enubias_bins = {-2.0, -0.6, -0.581, -0.5595, -0.538, -0.5165, -0.495, -0.4735, -0.452, -0.4305, -0.409, -0.3875, -0.366, -0.3445, -0.323, -0.3015, -0.28, -0.2585, -0.237, -0.2155, -0.194, -0.1725, -0.151, -0.1295, -0.108, -0.0865, -0.065, -0.0435, -0.022, 0.0, 0.1};
        
        static TH2D *tmplt_hist_observables = nullptr;
        static TH2D *tmplt_hist = nullptr;

        static TH1D *tmplt_1D_xsec = nullptr;
        static TH1D *tmplt_1D_observables = nullptr;
        std::vector<double> template1D_enu_bins ={0.0, 0.14492753623188406, 0.2898550724637681, 0.43478260869565216, 0.5797101449275363, 0.7246376811594203, 0.8695652173913044, 1.0144927536231884, 1.1594202898550723, 1.3043478260869565, 1.4492753623188405, 1.5942028985507247, 1.739130434782609, 1.8840579710144927, 2.028985507246377, 2.1739130434782608, 2.318840579710145, 2.463768115942029, 2.608695652173913, 2.753623188405797, 2.898550724637681, 3.0434782608695654, 3.188405797101449, 3.333333333333333, 3.4782608695652173, 3.6231884057971016, 3.7681159420289856, 3.9130434782608696, 4.057971014492754, 4.202898550724638, 4.347826086956522, 4.492753623188406, 4.63768115942029, 4.782608695652174, 4.927536231884058, 5.072463768115942, 5.217391304347826, 5.3623188405797105, 5.507246376811594, 5.652173913043478, 5.797101449275362, 5.9420289855072465, 6.086956521739131, 6.2318840579710145, 6.376811594202899, 6.521739130434783, 6.666666666666667, 6.811594202898551, 6.956521739130435, 7.101449275362319, 7.246376811594203, 7.391304347826087, 7.536231884057971, 7.681159420289856, 7.82608695652174, 7.971014492753624, 8.115942028985508, 8.260869565217392, 8.405797101449276, 8.55072463768116, 8.695652173913045, 8.840579710144928, 8.985507246376812, 9.130434782608696, 9.27536231884058, 9.420289855072464, 9.565217391304348, 9.710144927536232, 9.855072463768116, 10.0};
        std::vector<double> observables1D_enu_bins = {0.0, 0.20408163265306123, 0.40816326530612246, 0.6122448979591837, 0.8163265306122449, 1.0204081632653061, 1.2244897959183674, 1.4285714285714286, 1.6326530612244898, 1.836734693877551, 2.0408163265306123, 2.2448979591836737, 2.4489795918367347, 2.653061224489796, 2.857142857142857, 3.061224489795918, 3.2653061224489797, 3.469387755102041, 3.673469387755102, 3.8775510204081633, 4.081632653061225, 4.285714285714286, 4.4897959183673475, 4.693877551020408, 4.8979591836734695, 5.1020408163265305, 5.306122448979592, 5.510204081632653, 5.714285714285714, 5.918367346938776, 6.122448979591837, 6.326530612244898, 6.530612244897959, 6.7346938775510205, 6.938775510204081, 7.142857142857143, 7.346938775510204, 7.551020408163265, 7.755102040816327, 7.959183673469388, 8.16326530612245, 8.36734693877551, 8.571428571428571, 8.775510204081632, 8.979591836734694, 9.183673469387755, 9.387755102040817, 9.591836734693878, 9.795918367346939, 10.0};

        if(!tmplt_hist_observables){
          tmplt_hist_observables = new TH2D("tmplt_hist_obs", "", global_enu_bins.size()-1,global_enu_bins.data(), global_elep_bins.size()-1,global_elep_bins.data());
        }
        if(!tmplt_hist){
          tmplt_hist = new TH2D("tmplt_hist", "", tmplt_enu_bins.size()-1,tmplt_enu_bins.data(), tmplt_enubias_bins.size()-1,tmplt_enubias_bins.data());
        }
        if(!tmplt_1D_xsec){
          tmplt_1D_xsec = new TH1D("tmplt_hist_obs", "", template1D_enu_bins.size()-1,template1D_enu_bins.data());
        }
        if(!tmplt_1D_observables){
          tmplt_1D_observables = new TH1D("tmplt_hist_obs", "", observables1D_enu_bins.size()-1,observables1D_enu_bins.data());
        }
        
        //duneobj->template_global_bin_number[i] = tmplt_hist->FindFixBin(duneobj->rw_etru[i], duneobj->enu_bias[i]);
        duneobj->global_bin_number[i] = GetGenericBinningGlobalBinNumber(iSample, i); //--original version, will need this again
        
        //duneobj->global_bin_number[i] = tmplt_hist_observables->FindFixBin(duneobj->rw_erec[i], duneobj->rw_erec_lep[i]);         
        //duneobj->template_global_bin_number[i] = tmplt_hist->FindFixBin(duneobj->rw_etru[i], duneobj->enu_bias[i]);

        //duneobj->global_bin_number[i] = tmplt_1D_observables->FindFixBin(duneobj->rw_erec[i]);         
        //duneobj->template_global_bin_number[i] = tmplt_1D_xsec->FindFixBin(duneobj->rw_etru[i]);

        //std::cout << "global_bin_number " << duneobj->global_bin_number[i] << std::endl;
        //std::cout << "template_global_bin_number[i] " << duneobj->template_global_bin_number[i] << std::endl;

        //outfile << duneobj->global_bin_number[i] << " "<< duneobj->template_global_bin_number[i] << std::endl;

  }
  }

  std::cout << "Condition was satisfied " << conditionCounter << " times." << std::endl;

  
  _sampleFile->Close();
  outfile.close();
  return duneobj->nEvents;
}

double samplePDFDUNEBeamFD::ReturnKinematicParameter(int KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  double KinematicValue = -999;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = dunemcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = dunemcSamples[iSample].rw_erec[iEvent];
    break;
  case kRecoNeutrinoEnergy_shifted:
    KinematicValue = dunemcSamples[iSample].rw_erec_shifted[iEvent];
    break;
  case kTrueXPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_x[iEvent];
    break;
  case kTrueYPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_y[iEvent];
    break;
  case kTrueZPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_z[iEvent];
    break;
  case kCVNNumu:
    KinematicValue = dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
    break;
  case kCVNNue:
    KinematicValue = dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
    break;
  case kM3Mode:
    KinematicValue = dunemcSamples[iSample].mode[iEvent];
    break;
  case kOscChannel:
    KinematicValue = MCSamples[iSample].ChannelIndex;
    break;
  case kIsFHC:
    KinematicValue = isFHC;
    break;
  case kq0:
    KinematicValue = dunemcSamples[iSample].true_q0[iEvent];
    break;
  case kq3:
    KinematicValue = dunemcSamples[iSample].true_q3[iEvent];
    break;
  case k_pT:
    KinematicValue = dunemcSamples[iSample].rw_pt[iEvent];
    break;
  case k_pz:
    KinematicValue = dunemcSamples[iSample].rw_pz[iEvent];
    break;
  case k_global_bin_number:
    KinematicValue = dunemcSamples[iSample].global_bin_number[iEvent];
    break;
  case k_template_global_bin_number:
    KinematicValue = dunemcSamples[iSample].template_global_bin_number[iEvent];
    break;
  case kp_lep:
    KinematicValue = dunemcSamples[iSample].p_lep[iEvent];
    break;
  case ktheta_lep:
    KinematicValue = dunemcSamples[iSample].theta_lep[iEvent];
    break;
  case kELepRec:
    KinematicValue = dunemcSamples[iSample].rw_erec_lep[iEvent];
    break;
  case kEHadRec:
    KinematicValue = dunemcSamples[iSample].rw_erec_had[iEvent];
    break;
  case kERec_minus_Etrue:
    KinematicValue = dunemcSamples[iSample].enurec_minus_enutrue[iEvent];
    break;
  case kERecQE:
    KinematicValue = dunemcSamples[iSample].ERec_QE[iEvent];
    break;
  case kENuProxy_minus_Enutrue:
    KinematicValue = dunemcSamples[iSample].enu_proxy_minus_enutrue[iEvent];
    break;
  case kyRec:
    KinematicValue = dunemcSamples[iSample].rw_yrec[iEvent];
    break;
  case keHad_av:
    KinematicValue = dunemcSamples[iSample].eHad_av[iEvent];
    break;
  case kisCC:
    KinematicValue = dunemcSamples[iSample].rw_isCC[iEvent];
    break;
  case kisRelativeEnubias:
    KinematicValue = (dunemcSamples[iSample].relative_enu_bias[iEvent]);
    break;
  case kEnubias:
    KinematicValue = (dunemcSamples[iSample].enu_bias[iEvent]);
    break;
  default:
  MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", static_cast<int>(KinPar));
  MACH3LOG_ERROR("Was given a Kinematic Variable of {}", KinematicVariable);
  throw MaCh3Exception(__FILE__, __LINE__);
  }
  return KinematicValue;
}

double samplePDFDUNEBeamFD::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
 KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter)); 
 double KinematicValue = -999;
 
 switch(KinPar){
 case kTrueNeutrinoEnergy:
   KinematicValue = dunemcSamples[iSample].rw_etru[iEvent]; 
   break;
 case kRecoNeutrinoEnergy:
   KinematicValue = dunemcSamples[iSample].rw_erec[iEvent];
   break;
 case kRecoNeutrinoEnergy_shifted:
    KinematicValue = dunemcSamples[iSample].rw_erec_shifted[iEvent];
    break;
 case kTrueXPos:
   KinematicValue = dunemcSamples[iSample].rw_vtx_x[iEvent];
   break;
 case kTrueYPos:
   KinematicValue = dunemcSamples[iSample].rw_vtx_y[iEvent];
   break;
 case kTrueZPos:
   KinematicValue = dunemcSamples[iSample].rw_vtx_z[iEvent];
   break;
 case kCVNNumu:
   KinematicValue = dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
   break;
 case kCVNNue:
   KinematicValue = dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
   break;
 case kM3Mode:
    KinematicValue = dunemcSamples[iSample].mode[iEvent];
    break;
 case kOscChannel:
    KinematicValue = MCSamples[iSample].ChannelIndex;
    break;
 case kIsFHC:
    KinematicValue = isFHC;
    break;
 case kq0:
    KinematicValue = dunemcSamples[iSample].true_q0[iEvent];
    break;
 case kq3:
    KinematicValue = dunemcSamples[iSample].true_q3[iEvent];
    break;
 case k_pT:
    KinematicValue = dunemcSamples[iSample].rw_pt[iEvent];
    break;
 case k_pz:
    KinematicValue = dunemcSamples[iSample].rw_pz[iEvent];
    break;
 case k_global_bin_number:
    KinematicValue = dunemcSamples[iSample].global_bin_number[iEvent];
    break;
  case k_template_global_bin_number:
    KinematicValue = dunemcSamples[iSample].template_global_bin_number[iEvent];
    break;
 case kp_lep:
    KinematicValue = dunemcSamples[iSample].p_lep[iEvent];
    break;
 case ktheta_lep:
    KinematicValue = dunemcSamples[iSample].theta_lep[iEvent];
    break;
  case kELepRec:
    KinematicValue = dunemcSamples[iSample].rw_erec_lep[iEvent];
    break;
  case kEHadRec:
    KinematicValue = dunemcSamples[iSample].rw_erec_had[iEvent];
    break;
  case kERec_minus_Etrue:
    KinematicValue = dunemcSamples[iSample].enurec_minus_enutrue[iEvent];
    break;
  case kERecQE:
    KinematicValue = dunemcSamples[iSample].ERec_QE[iEvent];
    break;
  case kENuProxy_minus_Enutrue:
    KinematicValue = dunemcSamples[iSample].enu_proxy_minus_enutrue[iEvent];
    break;
  case kyRec:
    KinematicValue = dunemcSamples[iSample].rw_yrec[iEvent];
    break;
  case keHad_av:
    KinematicValue = dunemcSamples[iSample].eHad_av[iEvent];
    break;
  case kisCC:
    KinematicValue = dunemcSamples[iSample].rw_isCC[iEvent];
    break;
  case kisRelativeEnubias:
    KinematicValue = (dunemcSamples[iSample].relative_enu_bias[iEvent]);
    break;
  case kEnubias:
    KinematicValue = (dunemcSamples[iSample].enu_bias[iEvent]);
    break;
  default:
   MACH3LOG_ERROR("Did not recognise Kinematic Parameter type {}", KinematicParameter);
   throw MaCh3Exception(__FILE__, __LINE__);
 }
 
 return KinematicValue;
}


const double* samplePDFDUNEBeamFD::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
 KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter)); 
 double* KinematicValue = nullptr;
 
 switch(KinPar){
 case kTrueNeutrinoEnergy:
   KinematicValue = &dunemcSamples[iSample].rw_etru[iEvent]; 
   break;
 case kRecoNeutrinoEnergy:
   KinematicValue = &dunemcSamples[iSample].rw_erec[iEvent];
   break;
  case kRecoNeutrinoEnergy_shifted:
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
 case kq0:
   KinematicValue = &(dunemcSamples[iSample].true_q0[iEvent]);
   break;
 case kq3:
   KinematicValue = &(dunemcSamples[iSample].true_q3[iEvent]);
   break;
 case k_pT:
   KinematicValue = &(dunemcSamples[iSample].rw_pt[iEvent]);
   break;
 case k_pz:
   KinematicValue = &(dunemcSamples[iSample].rw_pz[iEvent]);
   break;
 case k_global_bin_number:
   KinematicValue = &(dunemcSamples[iSample].global_bin_number[iEvent]);
   break;
 case k_template_global_bin_number:
   KinematicValue = &(dunemcSamples[iSample].template_global_bin_number[iEvent]);
   break;
 case kp_lep:
   KinematicValue = &(dunemcSamples[iSample].p_lep[iEvent]);
   break;
 case ktheta_lep:
   KinematicValue = &(dunemcSamples[iSample].theta_lep[iEvent]);
   break;
 case kELepRec:
   KinematicValue = &(dunemcSamples[iSample].rw_erec_lep[iEvent]);
   break;
 case kEHadRec:
   KinematicValue = &(dunemcSamples[iSample].rw_erec_had[iEvent]);
   break;
 case kERec_minus_Etrue:
   KinematicValue = &(dunemcSamples[iSample].enurec_minus_enutrue[iEvent]);
   break;
 case kERecQE:
   KinematicValue = &(dunemcSamples[iSample].ERec_QE[iEvent]);
   break;
 case kENuProxy_minus_Enutrue:
    KinematicValue = &(dunemcSamples[iSample].enu_proxy_minus_enutrue[iEvent]);
    break;
 case kyRec:
    KinematicValue = &(dunemcSamples[iSample].rw_yrec[iEvent]);
    break;
 case keHad_av:
    KinematicValue = &(dunemcSamples[iSample].eHad_av[iEvent]);
    break;
  case kisCC:
    KinematicValue = &(dunemcSamples[iSample].rw_isCC[iEvent]);
    break;
  case kisRelativeEnubias:
    KinematicValue = &(dunemcSamples[iSample].relative_enu_bias[iEvent]);
    break;
  case kEnubias:
    KinematicValue = &(dunemcSamples[iSample].enu_bias[iEvent]);
    break;
 default:
   MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicParameter);
   throw MaCh3Exception(__FILE__, __LINE__);
 }
 
 return KinematicValue;
}

const double* samplePDFDUNEBeamFD::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  double* KinematicValue = nullptr;

  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &dunemcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &dunemcSamples[iSample].rw_erec[iEvent];
    break;
  case kRecoNeutrinoEnergy_shifted:
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
  case k_pz:
    KinematicValue = &(dunemcSamples[iSample].rw_pz[iEvent]);
  break;
  case k_pT:
    KinematicValue = &(dunemcSamples[iSample].rw_pt[iEvent]);
  break;
  case kq0:
    KinematicValue = &(dunemcSamples[iSample].true_q0[iEvent]);
  break;
  case kq3:
    KinematicValue = &(dunemcSamples[iSample].true_q3[iEvent]);
  break;
  case k_global_bin_number:
    KinematicValue = &(dunemcSamples[iSample].global_bin_number[iEvent]);
  break;
  case k_template_global_bin_number:
    KinematicValue = &(dunemcSamples[iSample].template_global_bin_number[iEvent]);
  break;
  case kp_lep:
    KinematicValue = &(dunemcSamples[iSample].p_lep[iEvent]);
  break;
  case ktheta_lep:
    KinematicValue = &(dunemcSamples[iSample].theta_lep[iEvent]);
  break;
  case kELepRec:
    KinematicValue = &(dunemcSamples[iSample].rw_erec_lep[iEvent]);
  break;
  case kEHadRec:
    KinematicValue = &(dunemcSamples[iSample].rw_erec_had[iEvent]);
  break;
  case kERec_minus_Etrue:
    KinematicValue = &(dunemcSamples[iSample].enurec_minus_enutrue[iEvent]);
  break;
  case kERecQE:
    KinematicValue = &(dunemcSamples[iSample].ERec_QE[iEvent]);
  break;
  case kENuProxy_minus_Enutrue:
    KinematicValue = &(dunemcSamples[iSample].enu_proxy_minus_enutrue[iEvent]);
  break;
  case kyRec:
    KinematicValue = &(dunemcSamples[iSample].rw_yrec[iEvent]);
    break;
  case keHad_av:
    KinematicValue = &(dunemcSamples[iSample].eHad_av[iEvent]);
    break;
  case kisCC:
    KinematicValue = &(dunemcSamples[iSample].rw_isCC[iEvent]);
    break;
  case kisRelativeEnubias:
    KinematicValue = &(dunemcSamples[iSample].relative_enu_bias[iEvent]);
    break;
  case kEnubias:
    KinematicValue = &(dunemcSamples[iSample].enu_bias[iEvent]);
    break;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", KinematicVariable);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

void samplePDFDUNEBeamFD::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);  

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent) {
     if (!(dunemcSamples[iSample].reco_numu[iEvent] &&
      (dunemcSamples[iSample].muon_contained[iEvent] || dunemcSamples[iSample].muon_tracker[iEvent]) &&
      dunemcSamples[iSample].Ehad_veto[iEvent] < 30)) {
  continue;  // Skip event if it fails the cut
    }
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
    case kRecoNeutrinoEnergy_shifted:
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

      //returnvec = AbiFuncc(config)
    case kTrueXPos:
    case kTrueYPos:
    case kTrueZPos:
    case kCVNNue:
    case kCVNNumu:
    case kq0:
    case kq3:
    case k_pT:
    case k_pz:
    case k_global_bin_number:
    case k_template_global_bin_number:
    case kp_lep:
    case ktheta_lep:
    case kELepRec:
      // ReturnVec.resize(XBinEdges.size());
      //   for (unsigned int bin_i=0;bin_i<XBinEdges.size();bin_i++) {ReturnVec[bin_i] = XBinEdges[bin_i];}
      //   break;
    case kEHadRec:
    case kERec_minus_Etrue:
    case kERecQE:
    case kENuProxy_minus_Enutrue:
    case kyRec:
    case keHad_av:
    case kisCC:
    case kisRelativeEnubias:
    case kEnubias:
      // ReturnVec.resize(XBinEdges.size());
      // for (unsigned int bin_i=0;bin_i<XBinEdges.size();bin_i++) {ReturnVec[bin_i] = XBinEdges[bin_i];}
      // break;
    default:
        MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", static_cast<int>(KinematicParameter));
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    return ReturnVec;

}

void samplePDFDUNEBeamFD::TotalEScaleND(const double * par,
                                        std::size_t iSample,
                                        std::size_t iEvent) {
  // MACH3LOG_INFO("iSample = %zu, sysActive = %d", 
  //             iSample, 
  //             ParHandler->IsParameterUsedInSample(kTotalEScaleND, iSample));

  //double shift = (*par) * dunemcSamples[iSample].rw_erec_had[iEvent];
  //dunemcSamples[iSample].rw_erec_shifted[iEvent] += shift;
// MACH3LOG_INFO("TotalEScaleND par = {}, shift = {}", *par, shift);
//  std::cout << "COUT TEST: par=" << *par
//           << " rw_erec_had=" << dunemcSamples[iSample].rw_erec_had[iEvent]
//           << " shift=" << ((*par) * dunemcSamples[iSample].rw_erec_had[iEvent])
//           << std::endl;
//   std::cout << "tot_escale_fd_pos = " << tot_escale_fd_pos << std::endl;
// dunemcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunemcSamples[iSample].rw_erec_had[iEvent];
// if (*par != 0)
//     std::cout << "Non-zero par: " << *par << std::endl;

if (dunemcSamples[iSample].rw_erec_shifted[iEvent] < 2.0 && *par != 0) {
    dunemcSamples[iSample].rw_erec_shifted[iEvent] = 4;
  }

}


void samplePDFDUNEBeamFD::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering functional parameters");
  // This function manually populates the map of functional parameters
  // Maps the name of the functional parameter to the pointer of the function

  // This is the part where we manually enter things
  // A lambda function has to be used so we can refer to a non-static member function
  RegisterIndividualFuncPar("TotalEScaleND",
                            kTotalEScaleND,
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleND(par, iSample, iEvent); });
  
  MACH3LOG_INFO("Finished registering functional parameters");
}


// void samplePDFDUNEBeamFD::resetShifts(int iEvent) {
//   // Reset the shifts to the original values
//   dunemcSamples[iSample].rw_erec_shifted[iEvent] = dunemcSamples[iSample].rw_erec[iEvent];
// }

