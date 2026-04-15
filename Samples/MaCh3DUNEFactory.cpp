#include "Samples/MaCh3DUNEFactory.h"


SampleHandlerBase* GetMaCh3DuneInstance(std::string SampleType, std::string SampleConfig, std::unique_ptr<ParameterHandlerGeneric>& xsec, const std::shared_ptr<OscillationHandler>&  BeamOscillator_, const std::shared_ptr<OscillationHandler>&  AtmOscillator_, BeamNDCov beamNDCov) {
  SampleHandlerBase *Sample;

  (void)beamNDCov;
  
  if (SampleType == "BeamFD") {
    Sample = new SampleHandlerBeamFD(SampleConfig, xsec.get(), BeamOscillator_);
  } else if (SampleType == "BeamND") {
    
    if (beamNDCov.NDCov_FHC == nullptr || beamNDCov.NDCov_RHC == nullptr || beamNDCov.NDCov_all == nullptr) {
      MACH3LOG_ERROR("NDCov objects are not defined");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    // TMatrixD* NDCov = nullptr;
    // manager* tempSampleManager = new manager(SampleConfig.c_str());
    // int isFHC = tempSampleManager->raw()["DUNESampleBools"]["isFHC"].as<int>();
    // if(isFHC) {NDCov = NDCov_FHC;}
    // else {NDCov = NDCov_RHC;}
    
    Sample = new SampleHandlerBeamND(SampleConfig, xsec.get(), beamNDCov); 
  } else if (SampleType == "Atm") {
    Sample = new SampleHandlerAtm(SampleConfig, xsec.get(), AtmOscillator_);
  } else if (SampleType == "BeamNDGAr") {
    Sample = new SampleHandlerBeamNDGAr(SampleConfig, xsec.get());
  }
  else {
    MACH3LOG_ERROR("Invalid SampleType: {} defined in {}", SampleType, SampleConfig);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return Sample;
}

std::vector<SampleHandlerBase*> MakeMaCh3DuneInstance(std::unique_ptr<Manager>& FitManager, std::unique_ptr<ParameterHandlerGeneric>&xsec){

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

  std::shared_ptr<OscillationHandler> AtmOscHandler;

  if (CheckNodeExists(FitManager->raw(), "General", "SharedNuOscillatorObject", "ATM") == false) {
    MACH3LOG_INFO("No SharedNuOscillatorObject for ATM specified, so not creating one");
  } else {
    std::string OscillatorConfig = Get<std::string>(FitManager->raw()["General"]["SharedNuOscillatorObject"]["ATM"], __FILE__, __LINE__);
    std::vector<const double*> OscParams = xsec->GetOscParsFromSampleName("ATM");
    AtmOscHandler = std::make_shared<OscillationHandler>(OscillatorConfig,true,OscParams,12);
  }

  std::shared_ptr<OscillationHandler> BeamOscHandler;

  if (CheckNodeExists(FitManager->raw(), "General", "SharedNuOscillatorObject", "Beam") == false) {
    MACH3LOG_INFO("No SharedNuOscillatorObject for Beam specified, so not creating one");
  } else {
    std::string OscillatorConfig = Get<std::string>(FitManager->raw()["General"]["SharedNuOscillatorObject"]["Beam"], __FILE__, __LINE__);
    std::vector<const double*> OscParams = xsec->GetOscParsFromSampleName("FD_");
    BeamOscHandler = std::make_shared<OscillationHandler>(OscillatorConfig,true,OscParams,12);
  }

  // ==========================================================
  
  BeamNDCov beamNDCov;

  // Get ND detector covariance matrix
  if (CheckNodeExists(FitManager->raw(), "General", "Systematics", "NDCovFile") ){
    std::string NDCovMatrixFile = FitManager->raw()["General"]["Systematics"]["NDCovFile"].as<std::string>();
    bool useCombinedNDCov = GetFromManager(FitManager->raw()["General"]["Systematics"]["UseCombinedNDCov"], true);

    TFile *NDCovFile = new TFile(NDCovMatrixFile.c_str(), "READ");

    TMatrixD* NDCov_FHC = NDCovFile->Get<TMatrixD>("nd_fhc_frac_cov");
    TMatrixD* NDCov_RHC = NDCovFile->Get<TMatrixD>("nd_rhc_frac_cov");
    TMatrixD* NDCov_all = NDCovFile->Get<TMatrixD>("nd_all_frac_cov");

    if (!(NDCov_FHC && NDCov_RHC && NDCov_all)) {
      MACH3LOG_ERROR("Could not find NDCov objects from file: {}",NDCovMatrixFile);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    beamNDCov.NDCov_FHC = NDCov_FHC;
    beamNDCov.NDCov_RHC = NDCov_RHC;
    beamNDCov.NDCov_all = NDCov_all;
    beamNDCov.useCombinedNDCov = useCombinedNDCov;

    NDCovFile->Close();
    delete NDCovFile;
  }

  //####################################################################################
  //Create SampleHandler Objs
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Loading DUNE samples..");
  std::vector<std::string> DUNESampleConfigs = FitManager->raw()["General"]["DUNESamples"].as<std::vector<std::string>>();

  for(unsigned int Sample_i = 0 ; Sample_i < DUNESampleConfigs.size() ; Sample_i++){

    Manager* tempSampleManager = new Manager(DUNESampleConfigs[Sample_i].c_str());
    std::string SampleType = tempSampleManager->raw()["SampleHandlerName"].as<std::string>();

    DUNEPdfs.push_back(GetMaCh3DuneInstance(SampleType, DUNESampleConfigs[Sample_i], xsec, BeamOscHandler, AtmOscHandler, beamNDCov));

    // Pure for debugging, lets us set which weights we don't want via the manager
#if DEBUG_DUNE_WEIGHTS==1
    DUNEPdfs.back()->setWeightSwitchOffVector(FitManager->getWeightSwitchOffVector());
    DUNEPdfs.back()->setXsecWeightSwitchOffVector(FitManager->getXsecWeightSwitchOffVector());
#endif
  }

  return DUNEPdfs;
}
