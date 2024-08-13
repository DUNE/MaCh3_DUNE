#include <TROOT.h>

#include "samplePDFDUNEBeamFDBase.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"

samplePDFDUNEBeamFDBase::samplePDFDUNEBeamFDBase(double pot, std::string mc_version, covarianceXsec* xsec_cov) : samplePDFBase(pot) { 
  std::cout << "- Using DUNE sample config in this file " << mc_version << std::endl;
  //ETA - safety feature so you can't pass a NULL xsec_cov
  if(xsec_cov == NULL){std::cerr << "[ERROR:] You've passed me a NULL xsec covariance matrix... I need this to setup splines!" << std::endl; throw;}
  init(pot, mc_version, xsec_cov);          
}

samplePDFDUNEBeamFDBase::~samplePDFDUNEBeamFDBase() {
}

void samplePDFDUNEBeamFDBase::init(double pot, std::string samplecfgfile, covarianceXsec *xsec_cov) {
  std::string mtupleprefix;
  std::string mtuplesuffix;
  std::string splineprefix;
  std::string splinesuffix;

  char* sample_char = (char*)samplecfgfile.c_str();
  //ETA - trying to read stuff from yaml file
  manager* SampleManager = new manager(sample_char);

  //Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  //Inputs
  mtupleprefix = SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  mtuplesuffix = SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  splineprefix = SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  splinesuffix = SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();

  //Binning
  BinningOpt = SampleManager->raw()["Binning"]["BinningOpt"].as<int>();  
  std::vector<double> sample_erec_bins = SampleManager->raw()["Binning"]["XVarBins"].as<std::vector<double>>();
  std::vector<double> sample_theta_bins = SampleManager->raw()["Binning"]["YVarBins"].as<std::vector<double>>();

  samplename = SampleManager->raw()["SampleName"].as<std::string>();

  std::vector<std::string> mtuple_files;
  std::vector<std::string> spline_files;
  std::vector<int> sample_vecno;
  std::vector<int> sample_oscnutype;
  std::vector<int> sample_nutype;
  std::vector<bool> sample_signal;
  
  //Loop over all the sub-samples
  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
	std::cout << "Found sub sample" << std::endl;
	mtuple_files.push_back(osc_channel["mtuplefile"].as<std::string>());
	spline_files.push_back(osc_channel["splinefile"].as<std::string>());
	sample_vecno.push_back(osc_channel["samplevecno"].as<int>());
	sample_nutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["nutype"].as<int>())));
	sample_oscnutype.push_back(PDGToProbs(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>())));
	sample_signal.push_back(osc_channel["signal"].as<bool>());
  }


  //Now loop over the kinematic cuts
  for ( auto const &SelectionCuts : SampleManager->raw()["SelectionCuts"]) {
	std::cout << "Looping over selection cuts " << std::endl;
	SelectionStr.push_back(SelectionCuts["KinematicStr"].as<std::string>());

	SelectionBounds.push_back(SelectionCuts["Bounds"].as<std::vector<double>>());
	std::cout << "Found cut on string " << SelectionCuts["KinematicStr"].as<std::string>() << std::endl;
	std::cout << "With bounds " << SelectionCuts["Bounds"].as<std::vector<double>>()[0] << " to " << SelectionCuts["Bounds"].as<std::vector<double>>()[1] << std::endl;
  }
  NSelections = SelectionStr.size();

  //Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  double erec_bin_edges[sample_erec_bins.size()];
  double theta_bin_edges[sample_theta_bins.size()];
  for(unsigned erec_i = 0 ; erec_i < sample_erec_bins.size() ; erec_i++){erec_bin_edges[erec_i] = sample_erec_bins[erec_i];}
  for(unsigned theta_i = 0 ; theta_i < sample_theta_bins.size() ; theta_i++){theta_bin_edges[theta_i] = sample_theta_bins[theta_i];}

  // create dunemc storage
  int nSamples = SampleManager->raw()["NSubSamples"].as<int>();
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }

  //Now down with yaml file for sample
  delete SampleManager;
  if(sample_oscnutype.size() != dunemcSamples.size()){std::cerr << "[ERROR:] samplePDFDUNEBeamFDBase::samplePDFDUNEBeamFDBase() - something went wrong either getting information from sample config" << std::endl; throw;}

  for(unsigned iSample=0 ; iSample < dunemcSamples.size() ; iSample++){
    setupDUNEMC((mtupleprefix+mtuple_files[iSample]+mtuplesuffix).c_str(), &dunemcSamples[sample_vecno[iSample]], pot, sample_nutype[iSample], sample_oscnutype[iSample], sample_signal[iSample]);
  }

  for (int i=0;i<nSamples;i++) {
    struct fdmc_base obj = fdmc_base();
    MCSamples.push_back(obj);
  }

  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    setupFDMC(&dunemcSamples[sample_vecno[iSample]], &MCSamples[sample_vecno[iSample]], (splineprefix+spline_files[iSample]+splinesuffix).c_str());
  }

  std::cout << "################" << std::endl;
  std::cout << "Setup FD MC   " << std::endl;
  std::cout << "################" << std::endl;

  // ETA - If xsec_cov hasn't been passed to the samplePDFDUNEBeamFDBase constructor then it's NULL
  // and the old funcitonality is kept
  // this calls this function in the core code
  // this needs to come after setupFDMC as otherwise MCSamples.splinefile will be NULL
  SetXsecCov(xsec_cov); 

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

  nFDDetectorSystPointers = funcParsIndex.size();
  FDDetectorSystPointers = std::vector<const double*>(nFDDetectorSystPointers);

  int func_it = 0;
  for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++func_it) {
	std::string name = funcParsNames.at(func_it);

	if (name == "TotalEScaleFD") {
	  tot_escale_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(tot_escale_fd_pos);
	}

	else if (name == "TotalEScaleSqrtFD") {
	  tot_escale_sqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(tot_escale_sqrt_fd_pos);
	}

	else if (name == "TotalEScaleInvSqrtFD") {
	  tot_escale_invsqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(tot_escale_invsqrt_fd_pos);
	}

	else if (name == "HadEScaleFD") {
	  had_escale_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(had_escale_fd_pos);
	}

	else if (name == "HadEScaleSqrtFD") {
	  had_escale_sqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(had_escale_sqrt_fd_pos);
	}

	else if (name == "HadEScaleInvSqrtFD") {
	  had_escale_invsqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(had_escale_invsqrt_fd_pos);
	}

	else if (name == "MuEScaleFD") {
	  mu_escale_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(mu_escale_fd_pos);
	}

	else if (name == "MuEScaleSqrtFD") {
	  mu_escale_sqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(mu_escale_sqrt_fd_pos);
	}

	else if (name == "MuEScaleInvSqrtFD") {
	  mu_escale_invsqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(mu_escale_invsqrt_fd_pos);
	}

	else if (name == "NEScaleFD") {
	  n_escale_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(n_escale_fd_pos);
	}

	else if (name == "NEScaleSqrtFD") {
	  n_escale_sqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(n_escale_sqrt_fd_pos);
	}

	else if (name == "NEScaleInvSqrtFD") {
	  n_escale_invsqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(n_escale_invsqrt_fd_pos);
	}

	else if (name == "EMEScaleFD") {
	  em_escale_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(em_escale_fd_pos);
	}

	else if (name == "EMEScaleSqrtFD") {
	  em_escale_sqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(em_escale_sqrt_fd_pos);
	}

	else if (name == "EMEScaleInvSqrtFD") {
	  em_escale_invsqrt_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(em_escale_invsqrt_fd_pos);
	}

	else if (name == "HadResFD") {
	  had_res_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(had_res_fd_pos);
	}

	else if (name == "MuResFD") {
	  mu_res_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(mu_res_fd_pos);
	}

	else if (name == "NResFD") {
	  n_res_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(n_res_fd_pos);
	}

	else if (name == "EMResFD") {
	  em_res_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(em_res_fd_pos);
	}
	else if (name == "CVNNumuFD") {
	  cvn_numu_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(cvn_numu_fd_pos);
	}
	else if (name == "CVNNueFD") {
	  cvn_nue_fd_pos = *it;
	  FDDetectorSystPointers[func_it] = xsec_cov->retPointer(cvn_nue_fd_pos);
	}

	else { 
	  std::cerr << "Found a functional parameter which wasn't specified in the xml | samplePDFDUNEBeamFDBase:" << name << std::endl;
	  throw;
	}
  }

  std::vector<std::string> spline_filepaths;
  std::cout << "Now setting up Splines" << std::endl;

  for(unsigned iSample=0 ; iSample < MCSamples.size() ; iSample++){
    spline_filepaths.push_back(splineprefix+spline_files[iSample]+splinesuffix);
  }

  splineFile = new splinesDUNE(xsec_cov);

    //////////////////////////////////
  // Now add samples to spline monolith
  //////////////////////////////////

  //ETA - do we need to do this here?
  //Can't we have one splineFile object for all samples as Dan does in atmospherics fit?
  //Then just add spline files to monolith?
  std::cout<<"Adding samples to spline monolith"<<std::endl;
  std::cout << "samplename is " << samplename << std::endl;
  std::cout << "BinningOpt is " << BinningOpt << std::endl;
  std::cout << "SampleDetID is " << SampleDetID << std::endl;
  std::cout << "spline_filepaths is of size " << spline_filepaths.size() << std::endl;

  splineFile->AddSample(samplename, BinningOpt, SampleDetID, spline_filepaths);

  // Print statements for debugging
  splineFile->PrintArrayDimension();
  splineFile->CountNumberOfLoadedSplines(false, 1);
  splineFile->TransferToMonolith();
  std::cout << "--------------------------------" <<std::endl;

  std::cout << "################" << std::endl;
  std::cout << "Setup FD splines   " << std::endl;
  std::cout << "################" << std::endl;

  SetupNormParameters();
  SetupWeightPointers();

  fillSplineBins();
  
  _sampleFile->Close();

  std::cout << "-------------------------------------------------------------------" <<std::endl;

  //The binning here is arbitrary, now we get info from cfg so the 
  //set1DBinning and set2Dbinning calls below will make the binning
  //to be what we actually want
  _hPDF1D = new TH1D("hErec_nue", "Reconstructed Energy", 200, 0 , 50.0);
  dathist = new TH1D("dat_nue","",200,0, 50.0); 
  _hPDF2D = new TH2D("blah","blah",15,0,50.0*1000,15,0,150);
  dathist2d = new TH2D("dat2d_nue","",15,0,1500,15,0,150);

  //ETA Don't forget the -1 on the size here, as it's number of bins not bin edges
  set1DBinning(sample_erec_bins.size()-1, erec_bin_edges);
  set2DBinning(sample_erec_bins.size()-1, erec_bin_edges, sample_theta_bins.size()-1, theta_bin_edges); 
}


void samplePDFDUNEBeamFDBase::setupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
	for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
	  //DB Setting total weight pointers
	  MCSamples[i].ntotal_weight_pointers[j] = 6;
	  MCSamples[i].total_weight_pointers[j] = new double*[MCSamples[i].ntotal_weight_pointers[j]];
	  MCSamples[i].total_weight_pointers[j][0] = &(dunemcSamples[i].pot_s);
	  MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].norm_s);
	  MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].osc_w[j]);
	  MCSamples[i].total_weight_pointers[j][3] = &(dunemcSamples[i].rw_berpaacvwgt[j]);
	  MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].flux_w[j]);
	  MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
	}
  }
}


void samplePDFDUNEBeamFDBase::setupDUNEMC(const char *sampleFile, dunemc_base *duneobj, double pot, int nutype, int oscnutype, bool signal, bool hasfloats) {
  
  // set up splines
  std::cout << "-------------------------------------------------------------------" << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;
  
  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree*)_sampleFile->Get("caf");

  if(_data){
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
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
  duneobj->pot_s = (pot)/norm->GetBinContent(2);
  std::cout<< "pot_s = " << duneobj->pot_s << std::endl;
  std::cout<< "norm_s = " << duneobj->norm_s << std::endl;
  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;
 
  std::cout << "signal: " << duneobj->signal << std::endl;
  std::cout << "nevents: " << duneobj->nEvents << std::endl;

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
  duneobj->xsec_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];

  duneobj->energyscale_w = new double[duneobj->nEvents];
  duneobj->mode = new int[duneobj->nEvents];
  duneobj->rw_lower_erec_1d = new double[duneobj->nEvents]; //lower erec bound for bin
  duneobj->rw_upper_erec_1d = new double[duneobj->nEvents]; //upper erec bound for bin
  duneobj->rw_lower_erec_2d = new double[duneobj->nEvents]; //lower erec bound for bin
  duneobj->rw_upper_erec_2d = new double[duneobj->nEvents]; //upper erec bound for bin

  // change so only points to one
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);

  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) // Loop through tree
    {
      _data->GetEntry(i);
      duneobj->rw_cvnnumu[i] = (double)_cvnnumu;
      duneobj->rw_cvnnue[i] = (double)_cvnnue;
      duneobj->rw_cvnnumu_shifted[i] = (double)_cvnnumu; 
      duneobj->rw_cvnnue_shifted[i] = (double)_cvnnue;
      if (iselike) {
	    duneobj->rw_erec[i] = (double)_erec_nue;
	    duneobj->rw_erec_shifted[i] = (double)_erec_nue; 
	    duneobj->rw_erec_had[i] = (double)_erec_had_nue;
	    duneobj->rw_erec_lep[i] = (double)_erec_lep_nue;
	  }
      else {
	    duneobj->rw_erec[i] = (double)_erec; 
	    duneobj->rw_erec_shifted[i] = (double)_erec; 
	    duneobj->rw_erec_had[i] = (double)_erec_had; 
	    duneobj->rw_erec_lep[i] = (double)_erec_lep; 
      }

	  duneobj->rw_eRecoP[i] = (double)_eRecoP; 
	  duneobj->rw_eRecoPip[i] = (double)_eRecoPip; 
	  duneobj->rw_eRecoPim[i] = (double)_eRecoPim; 
	  duneobj->rw_eRecoPi0[i] = (double)_eRecoPi0; 
	  duneobj->rw_eRecoN[i] = (double)_eRecoN; 

	  duneobj->rw_LepE[i] = (double)_LepE; 
	  duneobj->rw_eP[i] = (double)_eP; 
	  duneobj->rw_ePip[i] = (double)_ePip; 
	  duneobj->rw_ePim[i] = (double)_ePim; 
	  duneobj->rw_ePi0[i] = (double)_ePi0; 
	  duneobj->rw_eN[i] = (double)_eN; 

      duneobj->rw_etru[i] = (double)_ev;
      duneobj->rw_theta[i] = (double)_LepNuAngle;
      duneobj->rw_isCC[i] = _isCC;
      duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
      duneobj->rw_nuPDG[i] = _nuPDG;
      duneobj->rw_berpaacvwgt[i] = (double)_BeRPA_cvwgt;
      duneobj->rw_vtx_x[i] = (double)_vtx_x;
      duneobj->rw_vtx_y[i] = (double)_vtx_y;
      duneobj->rw_vtx_z[i] = (double)_vtx_z;

      //Assume everything is on Argon for now....
      duneobj->Target[i] = 40;
      
      duneobj->xsec_w[i] = 1.0;

      int mode= TMath::Abs(_mode);       
      duneobj->mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
 
      duneobj->energyscale_w[i] = 1.0;
      
      duneobj->flux_w[i] = 1.0;
    }
  
  _sampleFile->Close();
  std::cout << "Sample set up OK" << std::endl;  
}

double samplePDFDUNEBeamFDBase::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {

 KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter)); 
 double KinematicValue = -999;

 switch(KinPar){
   case kTrueNeutrinoEnergy:
	 KinematicValue = dunemcSamples[iSample].rw_etru[iEvent]; 
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
   default:
	 std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__ << " Did not recognise Kinematic Parameter type..." << std::endl;
	 throw;
 }

  return KinematicValue;
}

void samplePDFDUNEBeamFDBase::setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj, const char *splineFile)  {

  fdobj->nEvents = duneobj->nEvents;
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;
  fdobj->x_var = new double*[fdobj->nEvents];
  fdobj->y_var = new double*[fdobj->nEvents];
  fdobj->rw_etru = new double*[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];    
  fdobj->NomXBin = new int[fdobj->nEvents];
  fdobj->NomYBin = new int[fdobj->nEvents];
  fdobj->XBin = new int [fdobj->nEvents];
  fdobj->YBin = new int [fdobj->nEvents];;	 
  fdobj->rw_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_lower_lower_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->rw_upper_upper_xbinedge = new double [fdobj->nEvents];
  fdobj->mode = new int*[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents]; 
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->nxsec_norm_pointers = new int[fdobj->nEvents];
  fdobj->xsec_norm_pointers = new const double**[fdobj->nEvents];
  fdobj->xsec_norms_bins = new std::list< int >[fdobj->nEvents];
  fdobj->xsec_w = new double[fdobj->nEvents];
  fdobj->flux_w = new double[fdobj->nEvents];
  fdobj->osc_w = new double[fdobj->nEvents];
  fdobj->isNC = new bool[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double**[fdobj->nEvents];
  fdobj->ntotal_weight_pointers = new int[fdobj->nEvents];
  fdobj->total_weight_pointers = new double**[fdobj->nEvents];
  fdobj->Target = new int*[fdobj->nEvents];

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
	fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
	fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
	fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
	//ETA - set these to a dummy value to begin with, these get set in set1DBinning or set2DBinning
	fdobj->NomXBin[iEvent] = -1;
	fdobj->NomYBin[iEvent] = -1;
	fdobj->XBin[iEvent] = -1;
	fdobj->YBin[iEvent] = -1;	 
	fdobj->rw_lower_xbinedge[iEvent] = -1;
	fdobj->rw_lower_lower_xbinedge[iEvent] = -1;
	fdobj->rw_upper_xbinedge[iEvent] = -1;
	fdobj->rw_upper_upper_xbinedge[iEvent] = -1;
	fdobj->xsec_w[iEvent] = 1.0;
	fdobj->osc_w[iEvent] = 1.0;
	fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
	fdobj->flux_w[iEvent] = duneobj->flux_w[iEvent];

	//ETA - this is where the variables that you want to bin your samples in are defined
	//If you want to bin in different variables this is where you put it for now
	switch(BinningOpt){
	  case 0:
	  case 1:
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec_shifted[iEvent]);
		fdobj->y_var[iEvent] = &(duneobj->dummy_y);//ETA - don't think we even need this as if we have a 1D sample we never need this, just not sure I like an unitialised variable in fdmc struct? 
		break;
	  case 2:
		//Just point to xvar to the address of the variable you want to bin in
		//This way we don't have to update both fdmc and skmc when we apply shifts
		//to variables we're binning in
		fdobj->x_var[iEvent] = &(duneobj->rw_erec_shifted[iEvent]);
		fdobj->y_var[iEvent] = &(duneobj->rw_theta[iEvent]);
		break;
	  default:
		std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__ << " unrecognised binning option" << BinningOpt << std::endl;
		throw;
		break;
	}

  }
}

void samplePDFDUNEBeamFDBase::applyShifts(int iSample, int iEvent) {

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
}

std::vector<double> samplePDFDUNEBeamFDBase::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  std::cout << "ReturnKinematicVarBinning" << std::endl;
  std::vector<double> binningVector;
  return binningVector;
}
