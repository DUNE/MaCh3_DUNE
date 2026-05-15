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

  //####################################################################################
  //Create sample handler + parameter_handler objects
  auto [param_handler, samples] = MaCh3DuneFactory(FitManager);

  //Some place to store the histograms
  std::vector<TH1*> PredictionHistograms;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();


  for (auto handler : samples) {
    for (unsigned iSample = 0; iSample < handler->GetNSamples(); ++iSample) {

      std::string name = handler->GetSampleTitle(iSample);
      sample_names.push_back(name);
      TString NameTString = TString(name.c_str());

      handler->Reweight();

      TH1* data_hist;
      bool load_data = GetFromManager(FitManager->raw()["General"]["LoadDataFromFile"], false);
      if(load_data){
      	std::vector<std::string> data_hist_vec = Get<std::vector<std::string>>(FitManager->raw()["General"]["DataHistName"], __FILE__, __LINE__);
        auto data_file_name = Get<std::string>(FitManager->raw()["General"]["DataFileName"], __FILE__, __LINE__);

      	std::string data_hist_name = data_hist_vec[iSample];

      	MACH3LOG_INFO("Loaded data: {} from {} for {}", data_hist_name, data_file_name, name);
        auto data_file = TFile::Open(data_file_name.c_str(), "OPEN");
        if(!data_file || data_file->IsZombie()){
          MACH3LOG_ERROR("Cannot find data file {}", data_file_name);
          throw MaCh3Exception(__FILE__, __LINE__);
        }

        data_hist = data_file->Get<TH1>(data_hist_name.c_str());
        if(!data_hist){
          MACH3LOG_ERROR("Cannot find data histogram: {} in {}", data_hist_name, data_file_name);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
      else{
        data_hist = static_cast<TH1*>(handler->GetMCHist(iSample)->Clone(NameTString+"_DataHist"));
      }

      PredictionHistograms.push_back(data_hist);

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

  //###########################################################################################################
  //MCMC

  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);
  //Start chain from random position unless continuing a chain
  if(!StartFromPreviousChain){
    if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false)) {
      param_handler->ThrowParameters();
    }
  }

  //Add systematic objects
  MaCh3Fitter->AddSystObj(param_handler.get());

  if (StartFromPreviousChain) {
    std::string PreviousChainPath = FitManager->raw()["General"]["PosFileName"].as<std::string>();
    MACH3LOG_INFO("MCMC getting starting position from: {}",PreviousChainPath);
    MaCh3Fitter->StartFromPreviousFit(PreviousChainPath);
  }

  //Add samples
  for(auto Sample : samples){
    MaCh3Fitter->AddSampleHandler(Sample);
  }


  //Run fit
  MaCh3Fitter->RunMCMC();

  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", M3::Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", M3::Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
