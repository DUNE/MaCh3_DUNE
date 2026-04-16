#include "Samples/MaCh3DUNEFactory.h"

// ###############################################################
SampleHandlerBase* GetMaCh3DuneInstance(std::string SampleType, std::string SampleConfig, std::unique_ptr<ParameterHandlerGeneric>& param_handler, const std::shared_ptr<OscillationHandler>&  BeamOscillator_, const std::shared_ptr<OscillationHandler>&  AtmOscillator_, BeamNDCov beamNDCov) {
// ###############################################################
  SampleHandlerBase *Sample;

  (void)beamNDCov;
  
  if (SampleType == "BeamFD") {
    Sample = new SampleHandlerBeamFD(SampleConfig, param_handler.get(), BeamOscillator_);
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
    
    Sample = new SampleHandlerBeamND(SampleConfig, param_handler.get(), beamNDCov); 
  } else if (SampleType == "Atm") {
    Sample = new SampleHandlerAtm(SampleConfig, param_handler.get(), AtmOscillator_);
  } else if (SampleType == "BeamNDGAr") {
    Sample = new SampleHandlerBeamNDGAr(SampleConfig, param_handler.get());
  }
  else {
    MACH3LOG_ERROR("Invalid SampleType: {} defined in {}", SampleType, SampleConfig);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return Sample;
}

std::shared_ptr<OscillationHandler> SetupOscillationHandler(std::unique_ptr<Manager> &FitManager, std::unique_ptr<ParameterHandlerGeneric>& param_handler,
                                                            const std::string& osc_config_file, const std::string& sample_name){

  if (!CheckNodeExists(FitManager->raw(), "General", "SharedNuOscillatorObject", osc_config_file))
  {
    MACH3LOG_INFO("No SharedNuOscillatorObject for {} found, so not creating one", osc_config_file);
    return nullptr;
  }

  auto osc_config = Get<std::string>(FitManager->raw()["General"]["SharedNuOscillatorObject"][osc_config_file], __FILE__, __LINE__);
  auto osc_params = param_handler->GetOscParsFromSampleName(sample_name);
  if (osc_params.size()==0){
    MACH3LOG_ERROR("Tried to set up oscillator {} but found 0 oscillation parameters!", sample_name);
    throw MaCh3Exception("Oscillator has no oscillation parameters" __FILE__, __LINE__);
  }

  return std::make_shared<OscillationHandler>(osc_config, true, osc_params, 12);
}

// ###############################################################
BeamNDCov SetupBeamNDCov(std::unique_ptr<Manager> &FitManager)
// ###############################################################
{

  BeamNDCov beam_nd_cov;

  if (!CheckNodeExists(FitManager->raw(), "General", "Systematics", "NDCovFile"))
  {
    return beam_nd_cov;
  }

  std::string NDCovMatrixFile = FitManager->raw()["General"]["Systematics"]["NDCovFile"].as<std::string>();
  bool useCombinedNDCov = GetFromManager(FitManager->raw()["General"]["Systematics"]["UseCombinedNDCov"], true);

  auto NDCovFile = M3::Open(NDCovMatrixFile, "READ", __FILE__, __LINE__);

  TMatrixD *NDCov_FHC = NDCovFile->Get<TMatrixD>("nd_fhc_frac_cov");
  TMatrixD *NDCov_RHC = NDCovFile->Get<TMatrixD>("nd_rhc_frac_cov");
  TMatrixD *NDCov_all = NDCovFile->Get<TMatrixD>("nd_all_frac_cov");

  if (!(NDCov_FHC && NDCov_RHC && NDCov_all))
  {
    MACH3LOG_ERROR("Could not find NDCov objects from file: {}", NDCovMatrixFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  beam_nd_cov.NDCov_FHC = NDCov_FHC;
  beam_nd_cov.NDCov_RHC = NDCov_RHC;
  beam_nd_cov.NDCov_all = NDCov_all;
  beam_nd_cov.useCombinedNDCov = useCombinedNDCov;

  NDCovFile->Close();
  return beam_nd_cov;
}

// ###############################################################
std::vector<SampleHandlerBase *> MaCh3DuneSampleFactory(std::unique_ptr<Manager> &FitManager, std::unique_ptr<ParameterHandlerGeneric> &param_handler)
// ###############################################################
{

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
  // Setup oscillation handlers
  auto AtmOscHandler = SetupOscillationHandler(FitManager, param_handler, "ATM", "ATM");
  auto BeamOscHandler = SetupOscillationHandler(FitManager, param_handler, "BeamFD",  "BeamFD_Sample");
  // ==========================================================

  auto beamNDCov = SetupBeamNDCov(FitManager);

    // ####################################################################################
  // Create SampleHandler Objs
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Loading DUNE samples..");
  std::vector<std::string> DUNESampleConfigs = FitManager->raw()["General"]["DUNESamples"].as<std::vector<std::string>>();

  std::vector<SampleHandlerBase *> DUNEPdfs(DUNESampleConfigs.size());

  for (unsigned int Sample_i = 0; Sample_i < DUNESampleConfigs.size(); Sample_i++)
  {

    Manager* tempSampleManager = new Manager(DUNESampleConfigs[Sample_i].c_str());
    std::string SampleType = tempSampleManager->raw()["SampleHandlerName"].as<std::string>();

    auto sample = GetMaCh3DuneInstance(SampleType, DUNESampleConfigs[Sample_i], param_handler, BeamOscHandler, AtmOscHandler, beamNDCov);
    
    #if DEBUG_DUNE_WEIGHTS==1
    // Pure for debugging, lets us set which weights we don't want via the manager
    sample->setWeightSwitchOffVector(FitManager->getWeightSwitchOffVector());
    sample->setXsecWeightSwitchOffVector(FitManager->getXsecWeightSwitchOffVector());
    #endif
    DUNEPdfs[Sample_i] = sample;
  }

  return DUNEPdfs;
}
