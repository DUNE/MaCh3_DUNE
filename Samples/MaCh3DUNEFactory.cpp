#include "Samples/MaCh3DUNEFactory.h"

#ifdef BUILD_NDGAR
#include "Samples/SampleHandlerBeamNDGAr.h"
#else
//#include "Samples/SampleHandlerBeamFD.h"
//#include "Samples/SampleHandlerBeamND.h"
#include "Samples/SampleHandlerAtm.h"
#endif

SampleHandlerFD* GetMaCh3DuneInstance(std::string SampleType, std::string SampleConfig, ParameterHandlerGeneric* &xsec, TMatrixD* NDCov_FHC, TMatrixD* NDCov_RHC) {
  SampleHandlerFD *Sample;

  (void)NDCov_FHC;
  (void)NDCov_RHC;
  
#ifdef BUILD_NDGAR
  (void)NDCov_FHC;
  (void)NDCov_RHC;
  if (SampleType == "BeamNDGAr") {
    Sample = new SampleHandlerBeamNDGAr(SampleConfig, xsec);
  } else {
    MACH3LOG_ERROR("Invalid SampleType: {} defined in {}", SampleType, SampleConfig);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
#else
  /*
  if (SampleType == "BeamFD") {
    Sample = new SampleHandlerBeamFD(SampleConfig, xsec);
  } else if (SampleType == "BeamND") {
    
    if (NDCov_FHC == nullptr || NDCov_RHC == nullptr) {
      MACH3LOG_ERROR("NDCov objects are not defined");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    TMatrixD* NDCov = nullptr;
    manager* tempSampleManager = new manager(SampleConfig.c_str());
    int isFHC = tempSampleManager->raw()["DUNESampleBools"]["isFHC"].as<int>();
    if(isFHC) {NDCov = NDCov_FHC;}
    else {NDCov = NDCov_RHC;}
    
    Sample = new SampleHandlerBeamND(SampleConfig, xsec, NDCov);
    } else*/
  if (SampleType == "Atm") {
    Sample = new SampleHandlerAtm(SampleConfig, xsec);
  } else {
    MACH3LOG_ERROR("Invalid SampleType: {} defined in {}", SampleType, SampleConfig);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
#endif
  
  return Sample;
}

void MakeMaCh3DuneInstance(manager *FitManager, std::vector<SampleHandlerFD*> &DUNEPdfs, ParameterHandlerGeneric *&xsec){

  // there's a check inside the manager class that does this; left here for demonstrative purposes
  if (FitManager == nullptr) {
    MACH3LOG_ERROR("Didn't find a good config in input configuration");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //Check that you have specified some DUNE samples
  if(!FitManager->raw()["General"]["DUNESamples"]){
    MACH3LOG_ERROR("You didn't specify any DUNESample Configs to create samples from. Please add General:DUNESamples to your config");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // ==========================================================
  
  // Get inputted systematic parameters covariance matrices
  std::vector<std::string> xsecCovMatrixFile;
  if (CheckNodeExists(FitManager->raw(), "General", "Systematics", "XsecCovFile") ){
    xsecCovMatrixFile = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();
  } else {
    MACH3LOG_ERROR("Require General:Systematics:XsecCovFile node in {}, please add this to the file!", FitManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Setup the covariance matrices
  if(xsec == nullptr){
    xsec = new ParameterHandlerGeneric(xsecCovMatrixFile, "xsec_cov");
  }
  else{
    MACH3LOG_INFO("covariance Xsec has already been created so I am not re-initialising the object"); 
  }

  MACH3LOG_INFO("cov xsec setup");
  MACH3LOG_INFO("------------------------------");

  //read flat prior, fixed paramas from the config file
  std::vector<std::string> XsecFixParams = GetFromManager<std::vector<std::string>>(FitManager->raw()["General"]["Systematics"]["XsecFix"], {});

  // Fixed xsec parameters loop
  if (XsecFixParams.size() == 1 && XsecFixParams.at(0) == "All") {
    for (int j = 0; j < xsec->GetNumParams(); j++) {
      xsec->ToggleFixParameter(j);
    }
  } else {
    for (unsigned int j = 0; j < XsecFixParams.size(); j++) {
      xsec->ToggleFixParameter(XsecFixParams.at(j));
    }
  }
  MACH3LOG_INFO("xsec parameters loop done");

  xsec->SetParameters();
  xsec->SetStepScale(FitManager->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());

  std::vector<double> oscpars = FitManager->raw()["General"]["OscillationParameters"].as<std::vector<double>>();
  xsec->SetGroupOnlyParameters("Osc", oscpars);
  
  // ==========================================================
  
  TMatrixD* NDCov_FHC = nullptr;
  TMatrixD* NDCov_RHC = nullptr;

  // Get ND detector covariance matrix
  if (CheckNodeExists(FitManager->raw(), "General", "Systematics", "NDCovFile") ){
    std::string NDCovMatrixFile = FitManager->raw()["General"]["Systematics"]["NDCovFile"].as<std::string>();

    TFile *NDCovFile = new TFile(NDCovMatrixFile.c_str(), "READ");

    NDCov_FHC = NDCovFile->Get<TMatrixD>("nd_fhc_frac_cov");
    NDCov_RHC = NDCovFile->Get<TMatrixD>("nd_rhc_frac_cov");

    if (!(NDCov_FHC && NDCov_RHC)) {
      MACH3LOG_ERROR("Could not find NDCov objects from file: {}",NDCovMatrixFile);
    }

    NDCovFile->Close();
    delete NDCovFile;
  }

  //####################################################################################
  //Create SampleHandler Objs
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Loading DUNE samples..");
  std::vector<std::string> DUNESampleConfigs = FitManager->raw()["General"]["DUNESamples"].as<std::vector<std::string>>();

  for(unsigned int Sample_i = 0 ; Sample_i < DUNESampleConfigs.size() ; Sample_i++){

    manager* tempSampleManager = new manager(DUNESampleConfigs[Sample_i].c_str());
    std::string SampleType = tempSampleManager->raw()["SampleType"].as<std::string>();

    DUNEPdfs.push_back(GetMaCh3DuneInstance(SampleType, DUNESampleConfigs[Sample_i], xsec, NDCov_FHC, NDCov_RHC));

    // Pure for debugging, lets us set which weights we don't want via the manager
#if DEBUG_DUNE_WEIGHTS==1
    DUNEPdfs.back()->setWeightSwitchOffVector(FitManager->getWeightSwitchOffVector());
    DUNEPdfs.back()->setXsecWeightSwitchOffVector(FitManager->getXsecWeightSwitchOffVector());
#endif
  }

  return;
}
