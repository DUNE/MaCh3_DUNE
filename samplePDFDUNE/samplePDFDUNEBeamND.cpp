#include "samplePDFDUNEBeamND.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "TError.h"

//Here nullptr is passed instead of OscCov to prevent oscillation calculations from being performed for the ND Samples
samplePDFDUNEBeamND::samplePDFDUNEBeamND(std::string mc_version_, covarianceXsec* xsec_cov_,  TMatrixD* nd_cov_, covarianceOsc* osc_cov_=nullptr) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {
  OscCov = nullptr;
  
  if(!nd_cov_){
    MACH3LOG_ERROR("You've passed me a nullptr to a ND covarince matrix... ");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  NDCovMatrix = nd_cov_;

  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

samplePDFDUNEBeamND::~samplePDFDUNEBeamND() {
}

void samplePDFDUNEBeamND::Init() {
  dunendmcSamples.resize(nSamples,dunemc_base());
  
  IsFHC = SampleManager->raw()["DUNESampleBools"]["isFHC"].as<double>();
  iselike = SampleManager->raw()["DUNESampleBools"]["iselike"].as<bool>();
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
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];      
      MCSamples[i].total_weight_pointers[j][3] = &(dunendmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

int samplePDFDUNEBeamND::setupExperimentMC(int iSample) {

  // Set up the branches

  caf::StandardRecord* sr = new caf::StandardRecord();
  dunemc_base* duneobj = &dunendmcSamples[iSample];
  std::string FileName = mc_files[iSample];
  MACH3LOG_INFO("Reading File: {}",FileName);
  TFile* File = TFile::Open(FileName.c_str());
  if (!File || File->IsZombie()) {
    MACH3LOG_ERROR("Did not find File: {}",FileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  TTree* Tree = File->Get<TTree>("cafTree");
  if (!Tree){
    MACH3LOG_ERROR("Did not find Tree::cafTree in File: {}",FileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  Tree->SetBranchStatus("*", 1);
  Tree->SetBranchAddress("rec", &sr);

  if (IsFHC) { 
    duneobj->norm_s = (1e21/1.5e21);
  } else {
    duneobj->norm_s = (1e21/1.905e21);
  }
  duneobj->pot_s = (pot)/1e21;

  duneobj->nEvents = static_cast<int>(Tree->GetEntries());

  // Usual MaCh3 Variables
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
  // duneobj->rw_berpaacvwgt = new double[duneobj->nEvents];
  duneobj->nupdg = new int[duneobj->nEvents];
  duneobj->nupdgUnosc = new int[duneobj->nEvents];
  // duneobj->mode = new double[duneobj->nEvents];
  // duneobj->Target = new int[duneobj->nEvents];

  //Reco Particle Energy
  duneobj->rw_ndLAr_particle_eRecoMuon = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eRecoMuon->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_eRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eRecoP->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_eRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eRecoPip->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_eRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eRecoPim->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_eRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eRecoPi0->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_eRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eRecoN->reserve(7 * duneobj->nEvents);

  //True Particle Energy
  duneobj->rw_ndLAr_particle_eMuon = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eMuon->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_eP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eP->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ePip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ePip->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ePim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ePim->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ePi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ePi0->reserve(7 * duneobj->nEvents);
  duneobj->rw_ndLAr_particle_eN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_eN->reserve(7 * duneobj->nEvents);

  //Reco Particle Momentum
  duneobj->rw_ndLAr_particle_MomRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomRecoMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomRecoPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomRecoPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomRecoPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomRecoP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomRecoN->reserve(7*duneobj->nEvents);

  //True Particle Momentum
  duneobj->rw_ndLAr_particle_MomMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_MomN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_MomN->reserve(7*duneobj->nEvents);

  //Reco Particle Start Vertex
  duneobj->rw_ndLAr_particle_StartXvtxRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxRecoMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxRecoMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxRecoMu->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxRecoPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxRecoPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxRecoPip->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxRecoPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxRecoPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxRecoPim->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxRecoPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxRecoPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxRecoPi0->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxRecoP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxRecoP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxRecoP->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxRecoN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxRecoN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxRecoN->reserve(7*duneobj->nEvents);
  
  //True Particle Start vertex
  duneobj->rw_ndLAr_particle_StartXvtxMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxMu->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxPip->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxPim->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxPi0->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxP->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_StartXvtxN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartXvtxN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartYvtxN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartYvtxN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_StartZvtxN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_StartZvtxN->reserve(7*duneobj->nEvents);

  //Reco Particle End Vertex
  duneobj->rw_ndLAr_particle_EndXvtxRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxRecoMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxRecoMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxRecoMu->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxRecoPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxRecoPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxRecoPip->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxRecoPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxRecoPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxRecoPim->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxRecoPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxRecoPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxRecoPi0->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxRecoP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxRecoP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxRecoP->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxRecoN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxRecoN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxRecoN->reserve(7*duneobj->nEvents);

  //True Particle End Vertex
  duneobj->rw_ndLAr_particle_EndXvtxMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxMu->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxPip->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxPim->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxPi0->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxP->reserve(7*duneobj->nEvents);

  duneobj->rw_ndLAr_particle_EndXvtxN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndXvtxN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndYvtxN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndYvtxN->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_EndZvtxN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_EndZvtxN->reserve(7*duneobj->nEvents);

  // Reco Particle Angle
  duneobj->rw_ndLAr_particle_ThetaRecoMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaRecoMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaRecoPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaRecoPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaRecoPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaRecoPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaRecoPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaRecoPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaRecoP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaRecoP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaRecoN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaRecoN->reserve(7*duneobj->nEvents);

  // True Particle Angle
  duneobj->rw_ndLAr_particle_ThetaMu = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaMu->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaPip = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaPip->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaPim = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaPim->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaPi0 = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaPi0->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaP = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaP->reserve(7*duneobj->nEvents);
  duneobj->rw_ndLAr_particle_ThetaN = new std::vector<double>;
  duneobj->rw_ndLAr_particle_ThetaN->reserve(7*duneobj->nEvents);


  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    Tree->GetEntry(i);

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


/*
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    ND-LAr Particle Truth Information
    -----------------------------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------------------------
    */    
    
    int nprim = sr->mc.nu[0].nprim;
    for (int i = 0; i < nprim; i++) {
      int pdg = sr->mc.nu[0].prim[i].pdg;
      nparticlesinsample[iSample]++;
      double energy = static_cast<double>(sr->mc.nu[0].prim[i].p.E);
      double px = static_cast<double>(sr->mc.nu[0].prim[i].p.px);
      double py = static_cast<double>(sr->mc.nu[0].prim[i].p.py);
      double pz = static_cast<double>(sr->mc.nu[0].prim[i].p.pz);
      double start_pos_x = static_cast<double>(sr->mc.nu[0].prim[i].start_pos.x);
      double start_pos_y = static_cast<double>(sr->mc.nu[0].prim[i].start_pos.y);
      double start_pos_z = static_cast<double>(sr->mc.nu[0].prim[i].start_pos.z);
      double end_pos_x = static_cast<double>(sr->mc.nu[0].prim[i].end_pos.x);
      double end_pos_y = static_cast<double>(sr->mc.nu[0].prim[i].end_pos.y);
      double end_pos_z = static_cast<double>(sr->mc.nu[0].prim[i].end_pos.z);
      double momentum = sqrt(px * px + py * py + pz * pz);
      double theta = acos(pz / momentum);

      switch (pdg) {
      case 2212:
        duneobj->rw_ndLAr_particle_eP->push_back(energy);
        duneobj->rw_ndLAr_particle_MomRecoP->push_back(momentum);
        duneobj->rw_ndLAr_particle_ThetaRecoP->push_back(theta);
        duneobj->rw_ndLAr_particle_StartXvtxRecoP->push_back(start_pos_x);
        duneobj->rw_ndLAr_particle_StartYvtxRecoP->push_back(start_pos_y);
        duneobj->rw_ndLAr_particle_StartZvtxRecoP->push_back(start_pos_z);
        duneobj->rw_ndLAr_particle_EndXvtxRecoP->push_back(end_pos_x);
        duneobj->rw_ndLAr_particle_EndYvtxRecoP->push_back(end_pos_y);
        duneobj->rw_ndLAr_particle_EndZvtxRecoP->push_back(end_pos_z);
        break;
      case 211:
        duneobj->rw_ndLAr_particle_ePip->push_back(energy);
        duneobj->rw_ndLAr_particle_MomRecoPip->push_back(momentum);
        duneobj->rw_ndLAr_particle_ThetaRecoPip->push_back(theta);
        duneobj->rw_ndLAr_particle_StartXvtxRecoPip->push_back(start_pos_x);
        duneobj->rw_ndLAr_particle_StartYvtxRecoPip->push_back(start_pos_y);
        duneobj->rw_ndLAr_particle_StartZvtxRecoPip->push_back(start_pos_z);
        duneobj->rw_ndLAr_particle_EndXvtxRecoPip->push_back(end_pos_x);
        duneobj->rw_ndLAr_particle_EndYvtxRecoPip->push_back(end_pos_y);
        duneobj->rw_ndLAr_particle_EndZvtxRecoPip->push_back(end_pos_z);
        break;
      case -211:
        duneobj->rw_ndLAr_particle_ePim->push_back(energy);
        duneobj->rw_ndLAr_particle_MomRecoPim->push_back(momentum);
        duneobj->rw_ndLAr_particle_ThetaRecoPim->push_back(theta);
        duneobj->rw_ndLAr_particle_StartXvtxRecoPim->push_back(start_pos_x);
        duneobj->rw_ndLAr_particle_StartYvtxRecoPim->push_back(start_pos_y);
        duneobj->rw_ndLAr_particle_StartZvtxRecoPim->push_back(start_pos_z);
        duneobj->rw_ndLAr_particle_EndXvtxRecoPim->push_back(end_pos_x);
        duneobj->rw_ndLAr_particle_EndYvtxRecoPim->push_back(end_pos_y);
        duneobj->rw_ndLAr_particle_EndZvtxRecoPim->push_back(end_pos_z);
        break;
      case 111:
        duneobj->rw_ndLAr_particle_ePi0->push_back(energy);
        duneobj->rw_ndLAr_particle_MomRecoPi0->push_back(momentum);
        duneobj->rw_ndLAr_particle_ThetaRecoPi0->push_back(theta);
        duneobj->rw_ndLAr_particle_StartXvtxRecoPi0->push_back(start_pos_x);
        duneobj->rw_ndLAr_particle_StartYvtxRecoPi0->push_back(start_pos_y);
        duneobj->rw_ndLAr_particle_StartZvtxRecoPi0->push_back(start_pos_z);
        duneobj->rw_ndLAr_particle_EndXvtxRecoPi0->push_back(end_pos_x);
        duneobj->rw_ndLAr_particle_EndYvtxRecoPi0->push_back(end_pos_y);
        duneobj->rw_ndLAr_particle_EndZvtxRecoPi0->push_back(end_pos_z);
        break;
      case 2112:
        duneobj->rw_ndLAr_particle_eN->push_back(energy);
        duneobj->rw_ndLAr_particle_MomRecoN->push_back(momentum);
        duneobj->rw_ndLAr_particle_ThetaRecoN->push_back(theta);
        duneobj->rw_ndLAr_particle_StartXvtxRecoN->push_back(start_pos_x);
        duneobj->rw_ndLAr_particle_StartYvtxRecoN->push_back(start_pos_y);
        duneobj->rw_ndLAr_particle_StartZvtxRecoN->push_back(start_pos_z);
        duneobj->rw_ndLAr_particle_EndXvtxRecoN->push_back(end_pos_x);
        duneobj->rw_ndLAr_particle_EndYvtxRecoN->push_back(end_pos_y);
        duneobj->rw_ndLAr_particle_EndZvtxRecoN->push_back(end_pos_z);
        break;
      case 13:
        duneobj->rw_ndLAr_particle_eMuon->push_back(energy);
        duneobj->rw_ndLAr_particle_MomRecoMu->push_back(momentum);
        duneobj->rw_ndLAr_particle_ThetaRecoMu->push_back(theta);
        duneobj->rw_ndLAr_particle_StartXvtxRecoMu->push_back(start_pos_x);
        duneobj->rw_ndLAr_particle_StartYvtxRecoMu->push_back(start_pos_y);
        duneobj->rw_ndLAr_particle_StartZvtxRecoMu->push_back(start_pos_z);
        duneobj->rw_ndLAr_particle_EndXvtxRecoMu->push_back(end_pos_x);
        duneobj->rw_ndLAr_particle_EndYvtxRecoMu->push_back(end_pos_y);
        duneobj->rw_ndLAr_particle_EndZvtxRecoMu->push_back(end_pos_z);
        break;
      }

  }
  
  //_sampleFile->Close();
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

// Set the covariance matrix for this class
void samplePDFDUNEBeamND::setNDCovMatrix() {
  if (NDCovMatrix == NULL) {
    std::cerr << "Could not find ND Detector covariance matrix you provided to setCovMatrix" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  int nXBins = static_cast<int>(XBinEdges.size()-1);
  int nYBins = static_cast<int>(YBinEdges.size()-1);
  int covSize = nXBins*nYBins;

  if (covSize != NDCovMatrix->GetNrows()) {
    std::cerr << "Sample dimensions do not match ND Detector Covariance!" << std::endl;
    std::cerr << "Sample XBins * YBins = " << covSize  << std::endl;
    std::cerr << "ND Detector Covariance = " << NDCovMatrix->GetNrows() << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  
  std::vector<double> FlatCV;
  int iter = 0;

  // 2D -> 1D Array
  for (int xBin = 0; xBin < nXBins; xBin++) 
  {
    for (int yBin = 0; yBin < nYBins; yBin++) 
    {
        double CV = samplePDFFD_data[yBin][xBin];
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
double samplePDFDUNEBeamND::GetLikelihood()
{
  if (samplePDFFD_data == NULL) {
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

  int nXBins = static_cast<int>(XBinEdges.size()-1);
  int nYBins = static_cast<int>(YBinEdges.size()-1);

  int covSize = nXBins*nYBins;

  std::vector<double> FlatData;
  std::vector<double> FlatMCPred;

  //2D -> 1D 
  for (int xBin = 0; xBin < nXBins; xBin++) 
  {
    for (int yBin = 0; yBin < nYBins; yBin++) 
    {
        double MCPred = samplePDFFD_array[yBin][xBin];
        FlatMCPred.push_back(MCPred);

        double DataVal = samplePDFFD_data[yBin][xBin];
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
