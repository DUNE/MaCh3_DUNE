#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "Fitters//mcmc.h"
#include "Samples/MaCh3DUNEFactory.h"

int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //
  if(argc == 1){
    MACH3LOG_INFO("Usage: bin/Fit configs/config.yaml");
    return 1;
  }

  manager *FitManager = new manager(argv[1]);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  ParameterHandlerGeneric* xsec = nullptr;
  ParameterHandlerOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc); 

  //Some place to store the histograms
  std::vector<TH1*> PredictionHistograms;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  osc->SetParameters(FitManager->raw()["General"]["OscillationParameters"].as<std::vector<double>>());
  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {
    
    std::string name = DUNEPdfs[sample_i]->GetTitle();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    
    DUNEPdfs[sample_i] -> Reweight();
    if (DUNEPdfs[sample_i]->GetNDim() == 1){
      PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i] -> Get1DHist() -> Clone(NameTString+"_unosc")));
      DUNEPdfs[sample_i]->AddData(static_cast<TH1D*>(PredictionHistograms[sample_i]));
    }
      
    else if (DUNEPdfs[sample_i]->GetNDim() == 2){
      PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i] -> Get2DHist() -> Clone(NameTString+"_unosc")));
      DUNEPdfs[sample_i]->AddData(static_cast<TH2D*>(PredictionHistograms[sample_i]));
    }
    else {
      MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting"); 
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    
    
  }
  
  //Now print out some event rates, we'll make a nice latex table at some point 
  for (unsigned iPDF = 0; iPDF < DUNEPdfs.size() ; ++iPDF) {
    MACH3LOG_INFO("Integrals of nominal hists: ");
    MACH3LOG_INFO("{} : {}",sample_names[iPDF].c_str(),PredictionHistograms[iPDF]->Integral());
    MACH3LOG_INFO("--------------");
  }
  
  //###########################################################################################################
  //MCMC

  std::unique_ptr<mcmc> MaCh3Fitter = std::make_unique<mcmc>(FitManager);

  std::string throwmatrixfilename = GetFromManager<std::string>(FitManager->raw()["General"]["ThrowMatrixFile"], "");
  std::string throwmatrixname = GetFromManager<std::string>(FitManager->raw()["General"]["ThrowMatrixName"], "");
  if (throwmatrixfilename == "") {
    MACH3LOG_INFO("No throw matrix file specified, will throw from covariance matrix.");
  }
  else {
    TFile *throwmatrixfile = new TFile(throwmatrixfilename.c_str());
    if (throwmatrixfile->IsZombie()) {
      MACH3LOG_ERROR("Couldn't find {}", throwmatrixfilename);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    MACH3LOG_INFO("Found throw matrix file {}.", throwmatrixfilename);
    TMatrixDSym *throwmatrix = throwmatrixfile->Get<TMatrixDSym>(throwmatrixname.c_str());
    if (!throwmatrix) {
      MACH3LOG_ERROR("Couldn't find throw matrix {} in file {}", throwmatrixname, throwmatrixfilename);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    // Reset individual step scales to 1.0
    xsec->ResetIndivStepScale();
    // Keep global step scale flexible, but it should be 1
    double globalStepScale = FitManager->raw()["General"]["Systematics"]["XsecStepScale"].as<double>();
    if (globalStepScale != 1.0) {
      MACH3LOG_WARN("Global step scale is not 1.0, it is set to {}. This may cause issues when using adapted throw matrix.", globalStepScale);
    }
    xsec->SetStepScale(globalStepScale);
    MACH3LOG_WARN("I have set all the individual step scales to 1.0 since we are using an external throw matrix");
    xsec->SetThrowMatrix(throwmatrix);
    MACH3LOG_INFO("Set throw matrix from file {} with name {}",
                  throwmatrixfilename, throwmatrixname);
    // Print the throw matrix diagonals
    for (int i = 0; i < throwmatrix->GetNrows(); i++) {
      std::cout << (*throwmatrix)(i, i) << " ";
    }
    std::cout << std::endl;
  }

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);

  //Start chain from random position unless continuing a chain
  if(!StartFromPreviousChain){
    if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false)) {
      xsec->ThrowParameters();
    }
    osc->SetParameters();
    osc->ThrowParameters();
  }
  

  //Add systematic objects
  MaCh3Fitter->AddSystObj(osc);
  if (GetFromManager(FitManager->raw()["General"]["StatOnly"], false)){
    MACH3LOG_INFO("Running a stat-only fit so no systematics will be applied");
  }
  else {
    MaCh3Fitter->AddSystObj(xsec);
  }


  if (StartFromPreviousChain) {
    std::string PreviousChainPath = FitManager->raw()["General"]["PosFileName"].as<std::string>();
    MACH3LOG_INFO("MCMC getting starting position from: {}",PreviousChainPath);
    MaCh3Fitter->StartFromPreviousFit(PreviousChainPath);
  }
  
  //Add samples
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  
  //Run fit
  MaCh3Fitter->RunMCMC();

  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
