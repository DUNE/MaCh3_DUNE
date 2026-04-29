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

#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"

int main(int argc, char * argv[]) {

  auto FitManager = MaCh3ManagerFactory(argc, argv);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  ParameterHandlerGeneric* xsec = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  //Some place to store the histograms
  std::vector<TH1*> PredictionHistograms;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  // Adjust missing proton parameter only if it exists

    
  
    int idx = xsec->GetParIndex("MissingProtonFD");

    xsec->SetSingleParameter(idx, 0.2);
    xsec->SetFixParameter(idx);
      
  for (auto handler : DUNEPdfs) {
    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
    
      std::string name = handler->GetSampleTitle(iSample);
      sample_names.push_back(name);
      TString NameTString = TString(name.c_str());
      
      handler->Reweight();
      PredictionHistograms.push_back(static_cast<TH1*>(handler->GetMCHist(iSample)->Clone(NameTString+"_DataHist")));

      if (handler->GetNDim(iSample) == 1){
        handler->AddData(iSample, static_cast<TH1D*>(PredictionHistograms.back()));
      } else if (handler->GetNDim(iSample) == 2){
        handler->AddData(iSample, static_cast<TH2D*>(PredictionHistograms.back()));
      }
      
      else {
        MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting"); 
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }

      MACH3LOG_INFO("Integrals of nominal hists: ");
      MACH3LOG_INFO("{} : {}",name.c_str(),PredictionHistograms.back()->Integral());
      MACH3LOG_INFO("--------------");
    }
  }
 
    
    // xsec->ToggleFixParameter("MissingProtonFD");
    // xsec->SetPar(mpeIdx,0.2);
    
  //xsec->SetPar(xsec->GetParIndex("MissingProtonFD"), 0);
  //xsec->ToggleFixParameter("MissingProtonFD");
   for (unsigned iPDF = 0; iPDF < DUNEPdfs.size() ; ++iPDF) {
    MACH3LOG_INFO("Integrals of nominal hists: ");
    MACH3LOG_INFO("{} : {}",sample_names[iPDF].c_str(),PredictionHistograms[iPDF]->Integral());
    MACH3LOG_INFO("--------------");
  }
  
  //###########################################################################################################
  //MCMC

  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);
  //Start chain from random position unless continuing a chain
  if (!StartFromPreviousChain) {
    if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false)) {

      const int nPars = xsec->GetNParameters();
      std::vector<double> startPars(nPars);

      for (int i = 0; i < nPars; ++i) {
        startPars[i] = xsec->GetParInit(i);  // this is 1.0
      }

      xsec->SetParameters(startPars);
      MACH3LOG_INFO("Initialized xsec parameters to prior values");
    }
  }
 
    

  //Add systematic objects
  MaCh3Fitter->AddSystObj(xsec);
  

  if (StartFromPreviousChain) {
    std::string PreviousChainPath = FitManager->raw()["General"]["PosFileName"].as<std::string>();
    MACH3LOG_INFO("MCMC getting starting position from: {}",PreviousChainPath);
    MaCh3Fitter->StartFromPreviousFit(PreviousChainPath);
  }
  
  //Add samples
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

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
    xsec->ResetIndivStepScale();
    double globalStepScale = FitManager->raw()["General"]["Systematics"]["XsecStepScale"].as<double>();
    if (globalStepScale != 1.0) {
      MACH3LOG_WARN("Global step scale is not 1.0, it is set to {}. This may cause issues when using adapted throw matrix.", globalStepScale);
    }
    xsec->SetStepScale(globalStepScale);
    MACH3LOG_WARN("I have set all the individual step scales to 1.0 since we are using an external throw matrix");
    TMatrixTSym<double>* throwmatrix =
    dynamic_cast<TMatrixTSym<double>*>(
        throwmatrixfile->Get(throwmatrixname.c_str())
    );
    xsec->SetThrowMatrix(throwmatrix);
    MACH3LOG_INFO("Set throw matrix from file {} with name {}",
                  throwmatrixfilename, throwmatrixname);
    // Print the throw matrix diagonals
    for (int i = 0; i < throwmatrix->GetNrows(); i++) {
      std::cout << (*throwmatrix)(i, i) << " ";
    }
    std::cout << std::endl;
  }
   xsec->SetSingleParameter(idx, 0.0);
  xsec->SetFixParameter(idx);

  //Run fit
  MaCh3Fitter->RunMCMC();

  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
