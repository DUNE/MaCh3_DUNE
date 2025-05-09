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

#include "mcmc/mcmc.h"
#include "samplePDFDUNE/MaCh3DUNEFactory.h"

int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //
  if(argc == 1){
    MACH3LOG_INFO("Usage: bin/jointFitDUNEBeam configs/config.yaml");
    return 1;
  }

  manager *FitManager = new manager(argv[1]);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs

  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc); 

  //Some place to store the histograms
  std::vector<TH1*> PredictionHistograms;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {
    
    std::string name = DUNEPdfs[sample_i]->GetTitle();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    
    osc->setParameters();
    DUNEPdfs[sample_i] -> reweight();
    if (DUNEPdfs[sample_i]->GetNDim() == 1){
      PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc")));
      DUNEPdfs[sample_i]->addData(static_cast<TH1D*>(PredictionHistograms[sample_i]));
    }
      
    else if (DUNEPdfs[sample_i]->GetNDim() == 2){
      PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i] -> get2DHist() -> Clone(NameTString+"_unosc")));
      DUNEPdfs[sample_i]->addData(static_cast<TH2D*>(PredictionHistograms[sample_i]));
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

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);

  //Start chain from random position unless continuing a chain
  if(!StartFromPreviousChain){
    xsec->throwParameters();
    osc->throwParameters();
  }
  

  //Add systematic objects
  MaCh3Fitter->addSystObj(osc);
  if (GetFromManager(FitManager->raw()["General"]["StatOnly"], false)){
    MACH3LOG_INFO("Running a stat-only fit so no systematics will be applied");
  }
  else {
    MaCh3Fitter->addSystObj(xsec);
  }


  if (StartFromPreviousChain) {
    std::string PreviousChainPath = FitManager->raw()["General"]["PosFileName"].as<std::string>();
    MACH3LOG_INFO("MCMC getting starting position from: {}",PreviousChainPath);
    MaCh3Fitter->StartFromPreviousFit(PreviousChainPath);
  }
  
  //Add samples
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->addSamplePDF(Sample);
  }

  
  //Run fit
  MaCh3Fitter->runMCMC();

  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
