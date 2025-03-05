#include "samplePDFDUNE/MaCh3DUNEFactory.h"

#include "samplePDFDUNE/samplePDFDUNEBeamFD.h"
#include "samplePDFDUNE/samplePDFDUNEBeamND.h"
// #include "samplePDFDUNE/samplePDFDUNEBeamNDGar.h"
#include "samplePDFDUNE/samplePDFDUNEAtm.h"

samplePDFFDBase* GetMaCh3DuneInstance(std::string SampleType, std::string SampleConfig, covarianceXsec* &xsec, covarianceOsc* &osc) {

  samplePDFFDBase *FDSample;
  if (SampleType == "BeamFD") {
    FDSample = new samplePDFDUNEBeamFD(SampleConfig, xsec, osc);
  } else if (SampleType == "BeamND") {
    FDSample = new samplePDFDUNEBeamND(SampleConfig, xsec, osc);
  // } else if (SampleType == "BeamNDGar") {
  //   FDSample = new samplePDFDUNEBeamNDGar(SampleConfig, xsec);
  } else if (SampleType == "Atm") {
    FDSample = new samplePDFDUNEAtm(SampleConfig, xsec, osc);
  } else {
    MACH3LOG_ERROR("Invalid SampleType: {} defined in {}", SampleType, SampleConfig);
    throw MaCh3Exception(__FILE__, __LINE__);
  } 
  
  return FDSample;
}

void MakeMaCh3DuneInstance(manager *FitManager, std::vector<samplePDFFDBase*> &DUNEPdfs, covarianceXsec *&xsec, covarianceOsc *&osc){

  // there's a check inside the manager class that does this; left here for demonstrative purposes
  if (FitManager == nullptr) {
    MACH3LOG_ERROR("Didn't find a good config in input configuration");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  //Check that you have specified some DUNE samples
  if(!FitManager->raw()["General"]["DUNESamples"]){
    MACH3LOG_ERROR("You didn't specify any DUNESample configs to create samples from. Please add General:DUNESamples to your config");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
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
    xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
  }
  else{
    MACH3LOG_INFO("covariance Xsec has already been created so I am not re-initialising the object"); 
  }

  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  xsec->setParameters(XsecParVals);  
  xsec->setStepScale(FitManager->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
  
  MACH3LOG_INFO("cov xsec setup");
  MACH3LOG_INFO("------------------------------");

  // HH: Add a check to skip osc cov if OscCovFile is not specified
  std::vector<std::string> OscMatrixFile = GetFromManager<std::vector<std::string>>(FitManager->raw()["General"]["Systematics"]["OscCovFile"], {""});
  std::string  OscMatrixName = GetFromManager<std::string>(FitManager->raw()["General"]["Systematics"]["OscCovName"], "");
  std::vector<double> oscpars = GetFromManager<std::vector<double>>(FitManager->raw()["General"]["OscillationParameters"], {});
  std::string OscPars = "";
  bool useosc = oscpars.size() > 0;
  if (useosc){

    for (unsigned int i=0;i<oscpars.size();i++) {
      OscPars+=std::to_string(oscpars[i]);
      OscPars+=", ";
    }
    MACH3LOG_INFO("Oscillation Parameters being used: {} ", OscPars);
    
    osc = new covarianceOsc(OscMatrixFile,OscMatrixName.c_str());
    osc->setName("osc_cov");
    MACH3LOG_INFO("Osc cov setup");
    MACH3LOG_INFO("------------------------------");
  }  
  else {
    MACH3LOG_WARN("No OscillationParameters found, skipping oscillation parameters.");
  }

  // HH: moved this to be before fix or flat in accordance with the comments below:
  // Fill the parameter values with their nominal values
  // should _ALWAYS_ be done before overriding with fix or flat
  xsec->setParameters();
  if (useosc) {
    MACH3LOG_INFO("Setting oscillation parameters for osc");
    
    osc->setParameters(oscpars); 
  }
  
  // ==========================================================
  //read flat prior, fixed paramas from the config file
  std::vector<std::string> XsecFixParams   = GetFromManager<std::vector<std::string>>(FitManager->raw()["General"]["Systematics"]["XsecFix"], {""});
  
  // Fixed xsec parameters loop
  if (XsecFixParams.size() == 1 && XsecFixParams.at(0) == "All") {
    for (int j = 0; j < xsec->getNpars(); j++) {
      xsec->toggleFixParameter(j);
    }
  } else {
    for (unsigned int j = 0; j < XsecFixParams.size(); j++) {
      xsec->toggleFixParameter(XsecFixParams.at(j));
    }
  }
  MACH3LOG_INFO("xsec parameters loop done");
  
  
  //####################################################################################
  //Create samplePDFDUNE Objs
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Loading DUNE samples..");
  std::vector<std::string> DUNESampleConfigs = FitManager->raw()["General"]["DUNESamples"].as<std::vector<std::string>>();
  
  for(unsigned int Sample_i = 0 ; Sample_i < DUNESampleConfigs.size() ; Sample_i++){

    manager* tempSampleManager = new manager(DUNESampleConfigs[Sample_i].c_str());
    std::string SampleType = tempSampleManager->raw()["SampleType"].as<std::string>();
    
    DUNEPdfs.push_back(GetMaCh3DuneInstance(SampleType, DUNESampleConfigs[Sample_i], xsec, osc));

    // Pure for debugging, lets us set which weights we don't want via the manager
#if DEBUG_DUNE_WEIGHTS==1
    DUNEPdfs.back()->setWeightSwitchOffVector(FitManager->getWeightSwitchOffVector());
    DUNEPdfs.back()->setXsecWeightSwitchOffVector(FitManager->getXsecWeightSwitchOffVector());
#endif
  }

  // Adaptive MCMC stuff
  // HH: Removed the check for AdaptionOptions and left it to Fit.cpp or RunAdaptiveMCMC.cpp
  // if(FitManager->raw()["AdaptionOptions"])
  //   xsec->initialiseAdaption(FitManager->raw());
  
  return;
}
