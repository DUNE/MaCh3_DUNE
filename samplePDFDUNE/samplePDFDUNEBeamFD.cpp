#include "samplePDFDUNEBeamFD.h"
#include "TVector3.h"
#include <cstdlib>
#include <iostream>
#include <cassert>

#include <yaml-cpp/yaml.h>
#include <set>
#include <iostream>

samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  OscCov = nullptr;
  ExtractQ0Q3BinningFromYAML("/scratch/abipeake/MaCh3_DUNE/configs/CovObjs/q0q3_triangle_1275params.yaml");
  Extract1DBinningFromYAML("/scratch/abipeake/MaCh3_DUNE/configs/CovObjs/1D_Reconeutrinoenergy.yaml");
  SetupQ0Q3Hist();
  Setup1DHist();
  
  Initialise();
}

samplePDFDUNEBeamFD::~samplePDFDUNEBeamFD() {
}

TH2D* samplePDFDUNEBeamFD::GetQ0Q3FillingHist() const {
  return fQ0Q3Hist;
}
TH1D* samplePDFDUNEBeamFD::Get1DFillingHist() const {
  return f1DHist;
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

/*
void samplePDFDUNEBeamFD::ExtractQ0Q3BinningFromYAML(const std::string& yaml_file) {
  std::set<double> q0_set, q3_set;
  fxsec_BinDefs.clear();

  try {
    YAML::Node config = YAML::LoadFile(yaml_file);
    const auto& systs = config["Systematics"];
    int index = 0;

    for (const auto& systematic : systs) {
      const auto& sys = systematic["Systematic"];
      const auto& cuts = sys["KinematicCuts"];

      double q0_min = -1, q0_max = -1, q3_min = -1, q3_max = -1;

      for (const auto& cut : cuts) {
        if (cut["q0"]) {
          const auto& q0 = cut["q0"];
          if (q0.IsSequence() && q0.size() == 2) {
            q0_min = q0[0].as<double>();
            q0_max = q0[1].as<double>();
            q0_set.insert(q0_min);
            q0_set.insert(q0_max);
          }
        }
        if (cut["q3"]) {
          const auto& q3 = cut["q3"];
          if (q3.IsSequence() && q3.size() == 2) {
            q3_min = q3[0].as<double>();
            q3_max = q3[1].as<double>();
            q3_set.insert(q3_min);
            q3_set.insert(q3_max);
          }
        }
      }
      if (q0_max > q0_min && q3_max > q3_min) {
        fxsec_BinDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
      }

    }

    fQ0Edges.assign(q0_set.begin(), q0_set.end());
    fQ3Edges.assign(q3_set.begin(), q3_set.end());

    std::cout << "[samplePDFDUNE] Loaded " << fxsec_BinDefs.size() << " q0/q3 bins.\n";

  } catch (const std::exception& e) {
    std::cerr << "YAML parsing failed in samplePDFDUNE: " << e.what() << "\n";
  }
} */

std::pair<std::vector<double>, std::vector<double>> 
samplePDFDUNEBeamFD::ExtractQ0Q3BinningFromYAML(const std::string& yaml_file) {
  std::set<double> q0_set, q3_set;
  fxsec_BinDefs.clear();

  try {
    YAML::Node config = YAML::LoadFile(yaml_file);
    const auto& systs = config["Systematics"];
    int index = 0;

    for (const auto& systematic : systs) {
      const auto& sys = systematic["Systematic"];
      const auto& cuts = sys["KinematicCuts"];

      double q0_min = -1, q0_max = -1, q3_min = -1, q3_max = -1;

      for (const auto& cut : cuts) {
        if (cut["q0"]) {
          const auto& q0 = cut["q0"];
          if (q0.IsSequence() && q0.size() == 2) {
            q0_min = q0[0].as<double>();
            q0_max = q0[1].as<double>();
            q0_set.insert(q0_min);
            q0_set.insert(q0_max);
          }
        }
        if (cut["q3"]) {
          const auto& q3 = cut["q3"];
          if (q3.IsSequence() && q3.size() == 2) {
            q3_min = q3[0].as<double>();
            q3_max = q3[1].as<double>();
            q3_set.insert(q3_min);
            q3_set.insert(q3_max);
          }
        }
      }

      if (q0_max > q0_min && q3_max > q3_min) {
        fxsec_BinDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
      }
    }
    fQ0Edges.assign(q0_set.begin(), q0_set.end());
    fQ3Edges.assign(q3_set.begin(), q3_set.end());


    std::vector<double> q0_edges(q0_set.begin(), q0_set.end());
    std::vector<double> q3_edges(q3_set.begin(), q3_set.end());

    std::cout << "[samplePDFDUNE] Loaded " << fxsec_BinDefs.size() << " q0/q3 bins.\n";

    return {q0_edges, q3_edges};

  } catch (const std::exception& e) {
    std::cerr << "YAML parsing failed in samplePDFDUNE: " << e.what() << "\n";
    return {{}, {}};
  }
}


void samplePDFDUNEBeamFD::SetupQ0Q3Hist() {
  if (fQ0Q3Hist) delete fQ0Q3Hist;

  if (fQ0Edges.empty() || fQ3Edges.empty()) {
  std::cerr << "❌ Error: q0/q3 edges not set before calling SetupQ0Q3Hist()" << std::endl;
  throw std::runtime_error("Missing q0/q3 bin edges");
  }


  int nQ3Bins = static_cast<int>(fQ3Edges.size()) - 1;
  int nQ0Bins = static_cast<int>(fQ0Edges.size()) - 1;

  /*
  fQ0Q3Hist = new TH2D("fQ0Q3Hist", "q_{0} vs q_{3}; q_{3} [GeV]; q_{0} [GeV]",
                     nQ3Bins, &fQ3Edges[0],
                     nQ0Bins, &fQ0Edges[0]);*/
  
  fQ0Q3Hist = new TH2D("fQ0Q3Hist", "q_{0} vs q_{3}; q_{3} [GeV]; q_{0} [GeV]",
                     nQ3Bins, &fQ3Edges[0],
                     nQ0Bins, &fQ0Edges[0]);

  
  fQ0Q3Hist->SetDirectory(nullptr);
}

void samplePDFDUNEBeamFD::Setup1DHist() {
  if (f1DHist) delete f1DHist;
  if (f1DEdges.empty()) {
  std::cerr << " Error: q0 edges not set before calling Setup1DHist()" << std::endl;
  throw std::runtime_error("Missing q0 bin edges");
  }
  int n1DBins = static_cast<int>(f1DEdges.size()) - 1;
  f1DHist = new TH1D("fQ0Q3Hist", "q_{0} ; q_{3} [GeV]; q_{0} [GeV]",
                     n1DBins, &f1DEdges[0]);

  
  f1DHist->SetDirectory(nullptr);
}

/*
void samplePDFDUNEBeamFD::Extract1DBinningFromYAML(const std::string& yaml_file) {
  std::set<double> q0_set;
  fxsec_BinDefs_1D.clear();

  try {
    YAML::Node config = YAML::LoadFile(yaml_file);
    const auto& systs = config["Systematics"];
    int index = 0;

    for (const auto& systematic : systs) {
      const auto& sys = systematic["Systematic"];
      const auto& cuts = sys["KinematicCuts"];

      double q0_min = -1, q0_max = -1;

      for (const auto& cut : cuts) {
        if (cut["RecoNeutrinoEnergy"]) {
          const auto& q0 = cut["RecoNeutrinoEnergy"];
          if (q0.IsSequence() && q0.size() == 2) {
            q0_min = q0[0].as<double>();
            q0_max = q0[1].as<double>();
            q0_set.insert(q0_min);
            q0_set.insert(q0_max);
          }
        }
      }
      if (q0_max > q0_min) {
        fxsec_BinDefs_1D.push_back({index++, q0_min, q0_max});
      }
    }
    f1DEdges.assign(q0_set.begin(), q0_set.end());

    std::cout << "[samplePDFDUNE] Loaded " << fxsec_BinDefs.size() << " q0/q3 bins.\n";

  } catch (const std::exception& e) {
    std::cerr << "YAML parsing failed in samplePDFDUNE: " << e.what() << "\n";
  }
}*/

std::vector<double> samplePDFDUNEBeamFD::Extract1DBinningFromYAML(const std::string& yaml_file) {
  std::set<double> q0_set;
  fxsec_BinDefs_1D.clear();

  try {
    YAML::Node config = YAML::LoadFile(yaml_file);
    const auto& systs = config["Systematics"];
    int index = 0;

    for (const auto& systematic : systs) {
      const auto& sys = systematic["Systematic"];
      const auto& cuts = sys["KinematicCuts"];

      double q0_min = -1, q0_max = -1;

      for (const auto& cut : cuts) {
        if (cut["RecoNeutrinoEnergy"]) {
          const auto& q0 = cut["RecoNeutrinoEnergy"];
          if (q0.IsSequence() && q0.size() == 2) {
            q0_min = q0[0].as<double>();
            q0_max = q0[1].as<double>();
            q0_set.insert(q0_min);
            q0_set.insert(q0_max);
          }
        }
      }

      if (q0_max > q0_min) {
        fxsec_BinDefs_1D.push_back({index++, q0_min, q0_max});
      }
    }
    f1DEdges.assign(q0_set.begin(), q0_set.end());

    std::vector<double> q0_edges(q0_set.begin(), q0_set.end());

    std::cout << "[samplePDFDUNE] Loaded " << fxsec_BinDefs_1D.size() << " 1D q0 bins.\n";

    return q0_edges;

  } catch (const std::exception& e) {
    std::cerr << "YAML parsing failed in samplePDFDUNE: " << e.what() << "\n";
    return {};
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

   
  // Assume q0_edges and q3_edges were parsed already
  h_q0q3_filling = new TH2D("h_q0q3_filling", "q0 vs q3 from MC",
    static_cast<Int_t>(fQ3Edges.size() - 1), &fQ3Edges[0],
    static_cast<Int_t>(fQ0Edges.size() - 1), &fQ0Edges[0]);

  h_q0q3_filling->SetDirectory(nullptr);  // Prevent ROOT memory ownership issues

  // Assume q0_edges and q3_edges were parsed already
  h_1D_filling = new TH1D("h1D_filling", "q0 from MC",
    static_cast<Int_t>(f1DEdges.size() - 1), &f1DEdges[0]);

  h_1D_filling->SetDirectory(nullptr);  // Prevent ROOT memory ownership issues

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
  duneobj->p_lep = new double[duneobj->nEvents];
  duneobj->theta_lep = new double[duneobj->nEvents];
  duneobj->rw_yrec = new double[duneobj->nEvents];
  
  duneobj->enurec_minus_enutrue = new double[duneobj->nEvents];
  duneobj->ERec_QE = new double[duneobj->nEvents];
  duneobj->erec_proxy = new double[duneobj->nEvents];
  duneobj->erec_proxy_minus_enu = new double[duneobj->nEvents];
  duneobj->eHad_av = new double[duneobj->nEvents];///////////////////////new variable :) bit like the truth one in GENIE

  duneobj->reco_numu = new bool[duneobj->nEvents];
  duneobj->muon_contained = new bool[duneobj->nEvents];
  duneobj->muon_tracker = new bool[duneobj->nEvents];
  duneobj->Ehad_veto = new float[duneobj->nEvents];
  duneobj->isNC = new bool[duneobj->nEvents];

  _data->GetEntry(0);
  bool need_global_bin_numbers = (XVarStr == "global_bin_number");
  //FILL DUNE STRUCT
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

    duneobj->rw_erec[i] = (_erec);
    duneobj->rw_erec_shifted[i] = (_erec_nue); 
    duneobj->rw_yrec[i] = ((_erec -_Elep_reco)/_erec);
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
    //std::cout << "eHad_av = " << _eP + _ePip + _ePim + _ePi0 + _eOther + (_nipi0 *0.1349) << std::endl;
    

    duneobj->rw_erec_lep[i] = (_Elep_reco); 
    duneobj->rw_erec_had[i] = (_erec - _Elep_reco);
    duneobj->erec_proxy_minus_enu[i] = (_Elep_reco) + eHad_truth - _ev;

    //std::cout<< "_Elep_rec = " << (_Elep_reco) << std::endl;
    //std::cout<< "rw_erec_had[i] = " <<  (_erec -_Elep_reco) << std::endl;
    //std::cout<< "enu true[i] = " << ( _ev) << std::endl;
    //std::cout<< "yrec[i] = " << ((_erec -_Elep_reco)/_erec) << std::endl;
    //std::cout<< "erec_proxy_minus_enu[i] = " <<  (_Elep_reco) + eHad_truth - _ev << std::endl;
    duneobj->isNC[i] = _isNC;

    
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
      duneobj->global_bin_number[i] = GetGenericBinningGlobalBinNumber(iSample, i);
    } 

    //std::cout << "[DEBUG] Filling fQ0Q3Hist with _ev = " << _ev 
      //    << ", (_Elep_reco + eHad_truth - _ev) = " << (_Elep_reco + eHad_truth - _ev) << std::endl;

    //h_q0q3_filling->Fill((TVector3{_NuMomX, _NuMomY, _NuMomZ} -TVector3{_LepMomX, _LepMomY, _LepMomZ}).Mag() , (_ev - _LepE));
    fQ0Q3Hist->Fill( (_Elep_reco) + eHad_truth - _ev, _ev );
    f1DHist->Fill( _ev );
  //DB Grab the associated enum with the argument string
    //int bin = fQ0Q3Hist->FindBin((_Elep_reco + eHad_truth - _ev),_ev);
    //std::cout << "[DEBUG] Bin content after fill: " << fQ0Q3Hist->GetBinContent(bin) << std::endl;
    
    //std::cout << "[DEBUG] fQ3Edges: [" << fQ3Edges.front() << ", " << fQ3Edges.back() << "]\n";
    //std::cout << "[DEBUG] fQ0Edges: [" << fQ0Edges.front() << ", " << fQ0Edges.back() << "]\n";


    //std::cout << "Filling _eV=" << _ev << ", Enuproxy - Enu =" <<  (_Elep_reco) + eHad_truth - _ev << std::endl;

  }
  
  _sampleFile->Close();
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
  case kp_lep:
    KinematicValue = dunemcSamples[iSample].p_lep[iEvent];
    break;
  case ktheta_lep:
    KinematicValue = dunemcSamples[iSample].theta_lep[iEvent];
    break;
  case kELepRec:
    KinematicValue = dunemcSamples[iSample].rw_LepE[iEvent];
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
  case kERecProxy_minus_Enu:
    KinematicValue = dunemcSamples[iSample].erec_proxy_minus_enu[iEvent];
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
 case kp_lep:
    KinematicValue = dunemcSamples[iSample].p_lep[iEvent];
    break;
 case ktheta_lep:
    KinematicValue = dunemcSamples[iSample].theta_lep[iEvent];
    break;
  case kELepRec:
    KinematicValue = dunemcSamples[iSample].rw_LepE[iEvent];
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
  case kERecProxy_minus_Enu:
    KinematicValue = dunemcSamples[iSample].erec_proxy_minus_enu[iEvent];
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
 case kp_lep:
   KinematicValue = &(dunemcSamples[iSample].p_lep[iEvent]);
   break;
 case ktheta_lep:
   KinematicValue = &(dunemcSamples[iSample].theta_lep[iEvent]);
   break;
 case kELepRec:
   KinematicValue = &(dunemcSamples[iSample].rw_LepE[iEvent]);
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
 case kERecProxy_minus_Enu:
    KinematicValue = &(dunemcSamples[iSample].erec_proxy_minus_enu[iEvent]);
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
  case kp_lep:
    KinematicValue = &(dunemcSamples[iSample].p_lep[iEvent]);
  break;
  case ktheta_lep:
    KinematicValue = &(dunemcSamples[iSample].theta_lep[iEvent]);
  break;
  case kELepRec:
    KinematicValue = &(dunemcSamples[iSample].rw_LepE[iEvent]);
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
  case kERecProxy_minus_Enu:
    KinematicValue = &(dunemcSamples[iSample].erec_proxy_minus_enu[iEvent]);
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
 
void samplePDFDUNEBeamFD::applyShifts(int iSample, int iEvent) {

  (void)iSample;
  (void)iEvent;
   
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

  bool CCnumu {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14) && dunemcSamples[iSample].nupdgUnosc==2};
  bool CCnue {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==12) && dunemcSamples[iSample].nupdgUnosc==1};
  bool NotCCnumu {!(dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14)) && dunemcSamples[iSample].nupdgUnosc==2};


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
    case kp_lep:
    case ktheta_lep:
    case kELepRec:
    case kEHadRec:
    case kERec_minus_Etrue:
    case kERecQE:
    case kERecProxy_minus_Enu:
    case kyRec:
    case keHad_av:
    case kisCC:
      {
          std::vector<double> ReturnVec_1Dconfig = Extract1DBinningFromYAML("/scratch/abipeake/MaCh3_DUNE/configs/CovObjs/1D_Reconeutrinoenergy.yaml") ;
          std::pair<std::vector<double>, std::vector<double>>  ReturnVec_2Dconfig = ExtractQ0Q3BinningFromYAML("/scratch/abipeake/MaCh3_DUNE/configs/CovObjs//scratch/abipeake/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_10.0GeV_fullgrid_smallerbins.yaml");
          ReturnVec = ReturnVec_1Dconfig;
          break;
      }

    default:
      MACH3LOG_ERROR("Did not recognise Kinematic Parameter type: {}", static_cast<int>(KinematicParameter));
      throw MaCh3Exception(__FILE__, __LINE__);
  }      
  return ReturnVec;
}
