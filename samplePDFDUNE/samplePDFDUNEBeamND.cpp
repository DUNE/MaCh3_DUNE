#include <TROOT.h>

#include "samplePDFDUNEBeamND.h"
#include "StructsDUNE.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include <string>
#include "TMath.h"
#include "manager/manager.h"

samplePDFDUNEBeamND::samplePDFDUNEBeamND(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {  
  Initialise();
}

samplePDFDUNEBeamND::~samplePDFDUNEBeamND() {
}

void samplePDFDUNEBeamND::Init() {
  dunendmcSamples.resize(nSamples,dunemc_base());
  funcParsGrid.resize(nSamples);

  if (CheckNodeExists(SampleManager->raw(), "POT")) {
    pot = SampleManager->raw()["POT"].as<double>();
  } else{
    MACH3LOG_ERROR("POT not defined in {}, please add this!", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  // SetupFunctionalParameters();

  /*
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

  std::vector<std::string> funcParsNames = XsecCov->GetParsNamesFromDetID(SampleDetID, SystType::kFunc);
  std::vector<int> funcParsIndex = XsecCov->GetParsIndexFromDetID(SampleDetID, SystType::kFunc);

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
    SplineHandler = std::unique_ptr<splineFDBase>(new splinesDUNE(XsecCov));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

// === HH: Functional parameters ===
void samplePDFDUNEBeamND::TotalEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  // Debugging couts
  std::cout << "-------------------------------------------------------------------" <<std::endl;
  std::cout << "TotalEScaleND: " << *par << std::endl;
  std::cout << "duneendmcSamples size" << dunendmcSamples.size() << std::endl;
  std::cout << "iSample: " << iSample << std::endl;
  std::cout << "iEvent: " << iEvent << std::endl;
  std::cout << "rw_erec: " << dunendmcSamples[iSample].rw_erec[iEvent] << std::endl;
  std::cout << "rw_erec_shifted: " << dunendmcSamples[iSample].rw_erec_shifted[iEvent] << std::endl;
  std::cout << "rw_erec_had: " << dunendmcSamples[iSample].rw_erec_had[iEvent] << std::endl;
  std::cout << "rw_erec_lep: " << dunendmcSamples[iSample].rw_erec_lep[iEvent] << std::endl;
  std::cout << "rw_isCC: " << dunendmcSamples[iSample].rw_isCC[iEvent] << std::endl;
  std::cout << "rw_nuPDG: " << dunendmcSamples[iSample].rw_nuPDG[iEvent] << std::endl;
  std::cout << "iselike: " << iselike << std::endl;
  // Apply the scale
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_had[iEvent];
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
  if (!(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==14) && iselike) {
    dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_lep[iEvent];
  }
  std::cout << "New rw_erec_shifted: " << dunendmcSamples[iSample].rw_erec_shifted[iEvent] << std::endl;
  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void samplePDFDUNEBeamND::DebugShift(const double * par, std::size_t iSample, std::size_t iEvent) {
  if (dunendmcSamples[iSample].rw_erec[iEvent] < 2) {
  // Debugging couts
  std::cout << "-------------------------------------------------------------------" <<std::endl;
  std::cout << "Old rw_erec_shifted: " << dunendmcSamples[iSample].rw_erec_shifted[iEvent] << std::endl;
    dunendmcSamples[iSample].rw_erec_shifted[iEvent] = 4.0;
  std::cout << "New rw_erec_shifted: " << dunendmcSamples[iSample].rw_erec_shifted[iEvent] << std::endl;
  std::cout << "-------------------------------------------------------------------" <<std::endl;
  }
}

void samplePDFDUNEBeamND::RegisterFunctionalParameters() {
  std::cout << "Registering functional parameters" << std::endl;
  // This function manually populates the map of functional parameters
  // Maps the name of the functional parameter to the pointer of the function
  std::vector<std::string> funcParsNamesVec = {};
  
  // This is the part where we manually enter things
  funcParsNamesMap["TotalEScaleND"] = kTotalEScaleND;
  funcParsNamesVec.push_back("TotalEScaleND");
  // A lambda function has to be used so we can refer to a non-static member function
  funcParsFuncMap[kTotalEScaleND] = [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScale(par, iSample, iEvent); };

  funcParsNamesMap["DebugNothing"] = kDebugNothing;
  funcParsNamesVec.push_back("DebugNothing");
  funcParsFuncMap[kDebugNothing] = [this](const double * par, std::size_t iSample, std::size_t iEvent) {};

  funcParsNamesMap["DebugShift"] = kDebugShift;
  funcParsNamesVec.push_back("DebugShift");
  funcParsFuncMap[kDebugShift] = [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->DebugShift(par, iSample, iEvent); };

  // For every functional parameter in XsecCov that matches the name in funcParsNames, add it to the map
  for (std::vector<FuncPars>::iterator it = funcParsVec.begin(); it != funcParsVec.end(); ++it) {
    if (std::find(funcParsNamesVec.begin(), funcParsNamesVec.end(), (*it).name) != funcParsNamesVec.end()) {
      std::cout << "Adding functional parameter: " << (*it).name << std::endl;
      std::cout << "Adding it into funcParsMap with key: " << funcParsNamesMap[(*it).name] << std::endl;
      std::cout << "The address of the function is: " << &(*it) << std::endl;
      funcParsMap[funcParsNamesMap[(*it).name]] = &(*it);
    }
  }
  // funcParsMap[kTotalEScaleND]->name = "James";
  // std::cout << "I have changed the name at address: " << &(funcParsMap[kTotalEScaleND]->name) << " to " << funcParsMap[kTotalEScaleND]->name << std::endl;
}

void samplePDFDUNEBeamND::SetupFunctionalParameters() {
  std::cout << "Setting up functional parameters" << std::endl;
  funcParsVec = XsecCov->GetFuncParsFromDetID(SampleDetID);
  RegisterFunctionalParameters();
  // HH check: Not sure if SampleDetID is defined at this point
  // For each event, make a vector of pointers to the functional parameters
  for (std::size_t iSample = 0; iSample < dunendmcSamples.size(); ++iSample) {
    funcParsGrid[iSample].resize(static_cast<std::size_t>(dunendmcSamples[iSample].nEvents));
    for (std::size_t iEvent = 0; iEvent < static_cast<std::size_t>(dunendmcSamples[iSample].nEvents); ++iEvent) {
      // Now loop over the functional parameters and get a vector of enums corresponding to the functional parameters
      for (std::vector<FuncPars>::iterator it = funcParsVec.begin(); it != funcParsVec.end(); ++it) {
        // Check whether the interaction modes match
        bool ModeMatch = MatchCondition((*it).modes, static_cast<int>(std::round((dunendmcSamples[iSample].mode[iEvent]))));
        if (!ModeMatch) {
          MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent, (dunendmcSamples[iSample].mode[iEvent]), (*it).name);
          continue;
        }
        // Now check whether within kinematic bounds
        bool IsSelected = true;
        if ((*it).hasKinBounds) {
          for (std::size_t iKinPar = 0; iKinPar < (*it).KinematicVarStr.size(); ++iKinPar) {
            // Check lower bound
            if (ReturnKinematicParameter((*it).KinematicVarStr[iKinPar], iSample, iEvent) <= (*it).Selection[iKinPar][0]) {
              IsSelected = false;
              MACH3LOG_TRACE("Event {}, missed Kinematic var check ({}) for dial {}", iEvent, (*it).KinematicVarStr[iKinPar], (*it).name);
              continue;
            }
            // Check upper bound
            else if (ReturnKinematicParameter((*it).KinematicVarStr[iKinPar], iSample, iEvent) > (*it).Selection[iKinPar][1]) {
              MACH3LOG_TRACE("Event {}, missed Kinematic var check ({}) for dial {}", iEvent, (*it).KinematicVarStr[iKinPar], (*it).name);
              IsSelected = false;
              continue;
            }
          }
        }
        // Need to then break the event loop
        if(!IsSelected){
          MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}", iEvent, (*it).name);
          continue;
        }
        FuncParEnum funcparenum = funcParsNamesMap[(*it).name];
        // std::cout << "Adding functional parameter: " << (*it).name << " to funcParsGrid at iSample: " << iSample << " and iEvent: " << iEvent << std::endl;
        funcParsGrid.at(iSample).at(iEvent).push_back(funcparenum);
      }
    }
  }
  std::cout << "Finished setting up functional parameters" << std::endl;
}

// =================================

void samplePDFDUNEBeamND::SetupWeightPointers() {
  for (int i = 0; i < (int)dunendmcSamples.size(); ++i) {
    for (int j = 0; j < dunendmcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunendmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunendmcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunendmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

int samplePDFDUNEBeamND::setupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  int nupdgUnosc = sample_nupdgunosc[iSample];
  int nupdg = sample_nupdg[iSample];
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files.at(iSample));
  
  _sampleFile = TFile::Open(mc_files.at(iSample).c_str(), "READ");
  _data = (TTree*)_sampleFile->Get("caf");

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
  _data->SetBranchStatus("isFHC", 1);
  _data->SetBranchAddress("isFHC", &_isFHC);
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

  if (!IsRHC) { 
    duneobj->norm_s = (1e21/1.5e21);
  } else {
    duneobj->norm_s = (1e21/1.905e21);
  }
  duneobj->pot_s = (pot)/1e21;
  
  std::cout << "pot: " << pot << std::endl;
  std::cout << "pot_s: " << duneobj->pot_s << std::endl;
  std::cout << "norm_s: " << duneobj->norm_s << std::endl;

  duneobj->nEvents = _data->GetEntries();

  // HH: Downsampling by choosing the first X% of the events
  // HH TODO: Instead of choosing the first X% of the events, we should randomly sample X% of the events 
  if (CheckNodeExists(SampleManager->raw(), "Downsample")) {
    double downsample = SampleManager->raw()["Downsample"].as<double>();
    MACH3LOG_INFO("Downsample found in {}, will sample only {} of each file and scale POT correspondingly!", SampleManager->GetFileName(), downsample);
    duneobj->nEvents = (int)(duneobj->nEvents*downsample);
    duneobj->pot_s = duneobj->pot_s/downsample;
    MACH3LOG_INFO("New number of events: {}", duneobj->nEvents);
    MACH3LOG_INFO("New POT: {}", pot);
    MACH3LOG_INFO("New pot_s: {}", duneobj->pot_s);
    MACH3LOG_INFO("New norm_s: {}", duneobj->norm_s);
  } else{
    MACH3LOG_INFO("Downsample not defined in {}, continuing without downsampling!", SampleManager->GetFileName());
  }



  duneobj->rw_yrec = new double[duneobj->nEvents];
  duneobj->rw_erec_lep = new double[duneobj->nEvents];
  duneobj->rw_erec_had = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_isFHC = new double[duneobj->nEvents];
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
    
    duneobj->rw_erec[i] = (double)_erec;
    duneobj->rw_erec_shifted[i] = (double)_erec;
    duneobj->rw_erec_lep[i] = (double)_erec_lep;
    duneobj->rw_erec_had[i] = (double)(_erec - _erec_lep);
    duneobj->rw_yrec[i] = (double)((_erec - _erec_lep)/_erec);
    duneobj->rw_etru[i] = (double)_ev; // in GeV
    duneobj->rw_theta[i] = (double)_LepNuAngle;
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_isFHC[i] = (double)_isFHC;
    duneobj->rw_reco_q[i] = _reco_q;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = (double)_BeRPA_cvwgt;
    
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

    //Assume everything is on Argon for now....
    duneobj->Target[i] = 40;
    
    int mode= TMath::Abs(_mode);       
    duneobj->mode[i]=(double)GENIEMode_ToMaCh3Mode(mode, _isCC);
    
    duneobj->flux_w[i] = 1.0;
  }
  
  _sampleFile->Close();
  return duneobj->nEvents;
}


const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &dunendmcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kRecoQ:
    KinematicValue = &dunendmcSamples[iSample].rw_reco_q[iEvent];
    break;
  case kRecoNeutrinoEnergy:
    // HH: Changed this to match BeamFD
    KinematicValue = &dunendmcSamples[iSample].rw_erec_shifted[iEvent];
    break;
  case kIsFHC:
    KinematicValue = &dunendmcSamples[iSample].rw_isFHC[iEvent];
    break;
  case kRecoY:
    KinematicValue = &dunendmcSamples[iSample].rw_yrec[iEvent];
    break;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    std::cout << KinPar << ReturnStringFromKinematicParameter(KinPar) << std::endl;
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
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

  fdobj->SampleDetID = SampleDetID;
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
  for (std::vector<FuncParEnum>::iterator it = funcParsGrid.at(iSample).at(iEvent).begin(); it != funcParsGrid.at(iSample).at(iEvent).end(); ++it) {
    // Check if func exists
    if (funcParsMap.find(*it) == funcParsMap.end()) {
      MACH3LOG_ERROR("Functional parameter {} not found in map", *it);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    funcParsFuncMap[*it](XsecCov->retPointer(funcParsMap[*it]->index), iSample, iEvent);
  }
}


/*
void samplePDFDUNEBeamND::applyShifts(int iSample, int iEvent) {
  // reset erec back to original value
  
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] = dunendmcSamples[iSample].rw_erec[iEvent];

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

  bool CCnumu {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==NuPDG::kNumu && dunendmcSamples[iSample].rw_nuPDGunosc[iEvent]==NuPDG::kNumu};
  bool CCnue {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==NuPDG::kNue && dunendmcSamples[iSample].rw_nuPDGunosc[iEvent]==NuPDG::kNue};
  bool NotCCnumu {!(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==14) && dunendmcSamples[iSample].rw_nuPDGunosc[iEvent]==NuPDG::kNumu};

/*
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
}
*/

std::vector<double> samplePDFDUNEBeamND::ReturnKinematicParameterBinning(std::string KinematicParameterStr) 
{
  std::vector<double> binningVector;
  return binningVector;
}

int samplePDFDUNEBeamND::ReturnKinematicParameterFromString(std::string KinematicParameterStr) {
  if(KinematicParameterStr == "TrueNeutrinoEnergy") return kTrueNeutrinoEnergy;
  if(KinematicParameterStr == "RecoQ") return kRecoQ;
  if(KinematicParameterStr == "RecoNeutrinoEnergy") return kRecoNeutrinoEnergy;
  if(KinematicParameterStr == "IsFHC") return kIsFHC;
  if(KinematicParameterStr == "RecoY") return kRecoY;
  return -1;
}

std::string samplePDFDUNEBeamND::ReturnStringFromKinematicParameter(int KinematicParameter) {
  switch(KinematicParameter){
    case kTrueNeutrinoEnergy: return "TrueNeutrinoEnergy";
    case kRecoQ: return "RecoQ";
    case kRecoNeutrinoEnergy: return "RecoNeutrinoEnergy";
    case kIsFHC: return "IsFHC";
    case kRecoY: return "RecoY";
    default: return "";
  }
}
