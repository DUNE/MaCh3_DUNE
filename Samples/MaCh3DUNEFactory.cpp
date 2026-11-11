#include "Samples/MaCh3DUNEFactory.h"

// DUNE Handlers
#include "Samples/SampleHandlerBeamFD.h"
#include "Samples/SampleHandlerBeamND.h"
#include "Samples/SampleHandlerBeamNDGAr.h"
#include "Samples/SampleHandlerAtm.h"
#include "Samples/SampleHandlerPDSP.h"

// ###############################################################
SampleHandlerBase* GetMaCh3DuneInstance(std::string SampleType, std::string SampleConfig, std::unique_ptr<ParameterHandlerGeneric>& param_handler, const std::shared_ptr<OscillationHandler>&  BeamOscillator_, const std::shared_ptr<OscillationHandler>&  AtmOscillator_, BeamNDCov beamNDCov) {
// ###############################################################
  SampleHandlerBase *Sample;

  if (SampleType == "BeamFD") {
    Sample = new SampleHandlerBeamFD(SampleConfig, param_handler.get(), BeamOscillator_);
  } else if (SampleType == "BeamND") {

    if (beamNDCov.NDCov_FHC == nullptr || beamNDCov.NDCov_RHC == nullptr || beamNDCov.NDCov_all == nullptr) {
      MACH3LOG_ERROR("NDCov objects are not defined");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    Sample = new SampleHandlerBeamND(SampleConfig, param_handler.get(), beamNDCov);
  } else if (SampleType == "Atm") {
    Sample = new SampleHandlerAtm(SampleConfig, param_handler.get(), AtmOscillator_);
  } else if (SampleType == "BeamNDGAr") {
    Sample = new SampleHandlerBeamNDGAr(SampleConfig, param_handler.get());
  } else if (SampleType == "PDSP") {
    Sample = new SampleHandlerPDSP(SampleConfig, param_handler.get());
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
  std::vector<const M3::float_t*> par_values(osc_params.size());
  std::vector<std::string> nu_osc_name(osc_params.size());
  for(size_t ij = 0; ij < osc_params.size(); ij++){
    par_values[ij] = param_handler->RetPointer(osc_params[ij].index);
    nu_osc_name[ij] = osc_params[ij].NuOscName;
  }

  return std::make_shared<OscillationHandler>(osc_config, true, par_values, nu_osc_name, 12);
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
  bool useCombinedNDCov = GetFromManager(FitManager->raw()["General"]["Systematics"]["UseCombinedNDCov"], true, __FILE__, __LINE__);

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

bool UsingDataFromFile(std::unique_ptr<Manager>& FitManager) {
  const auto data_file_name = GetFromManager<std::string>(FitManager->raw()["General"]["DataFile"], "", __FILE__, __LINE__);
  return !data_file_name.empty();
}

std::pair<std::unique_ptr<ParameterHandlerGeneric>, std::vector<SampleHandlerBase*>> MaCh3DuneFactory(std::unique_ptr<Manager> &FitManager) {
  /// Generates a MaCh3 DUNE instance
  auto param_handler = MaCh3CovarianceFactory<ParameterHandlerGeneric>(FitManager.get(), "Xsec");

  if (CheckNodeExists(FitManager->raw(), "General", "OscillationParameters")){
    auto oscpars = Get<std::vector<double>>(FitManager->raw()["General"]["OscillationParameters"], __FILE__, __LINE__);
    param_handler->SetGroupOnlyParameters("Osc", oscpars);
  }

  auto samples = MaCh3DuneSampleFactory(FitManager, param_handler);

  if(samples.empty()){
    MACH3LOG_WARN("Cannot find any samples, doing prior-only fit");
  }

  if(!UsingDataFromFile(FitManager)){
    MACH3LOG_INFO("General.DataFile is not set, so the samples hold no data yet");
    return {std::move(param_handler), samples};
  }

  auto data_file_name = Get<std::string>(FitManager->raw()["General"]["DataFile"], __FILE__, __LINE__);

  auto data_file = M3::Open(data_file_name, "READ", __FILE__, __LINE__);
  if(!data_file){
    MACH3LOG_ERROR("Cannot open data file: {}", data_file_name);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for(auto& handler : samples){
    for (int iSample=0; iSample<handler->GetNSamples(); iSample++) {
      auto sample_name = handler->GetSampleName(iSample);
      // Histograms are the sample name with h at the front, e.g. hFD_FHC_numu_cc_1Recoil
      // Remove Beam from the sample name if it exists, as the histograms in the data file don't have this prefix
      if (sample_name.find("Beam") == 0) {
        sample_name = sample_name.substr(4);
      }
      std::string hist_name = "h" + sample_name;
      auto hist = data_file->Get<TH1>(hist_name.c_str());
      if(!hist){
        MACH3LOG_ERROR("Cannot find histogram {} in data file: {}", hist_name, data_file_name);
        throw MaCh3Exception(__FILE__, __LINE__); 
      }
      handler->AddData(iSample, hist);
      MACH3LOG_INFO("Data for {} taken from {}:{}, integral {:.2f}",
                    handler->GetSampleTitle(iSample), data_file_name, hist_name, hist->Integral());
    }

  }

  return {std::move(param_handler), samples};
}
