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

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

int main(int argc, char * argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  // 1D scan on by default, and 2D off
  const bool do_1d_llhscan = GetFromManager(FitManager->raw()["General"]["1DLLHScan"], true);
  const bool do_2d_llhscan = GetFromManager(FitManager->raw()["General"]["2DLLHScan"], false);

  if (!do_1d_llhscan && !do_2d_llhscan) {
    MACH3LOG_ERROR("Neither 1D or 2D llhscan enabled");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //###############################################################################################################################
  //Create sample handler + parameter_handler objects
  auto [param_handler, samples] = MaCh3DuneFactory(FitManager);

  //###############################################################################################################################
  //Perform reweight, print total integral, and set data

  std::vector<std::unique_ptr<TH1>> DUNEHists;
  for(auto handler : samples){
    for (unsigned iSample = 0; iSample < handler->GetNSamples(); ++iSample) {
      handler->Reweight();
      DUNEHists.push_back(M3::Clone(handler->GetMCHist(iSample)));
      MACH3LOG_INFO("Event rate for {} : {:<5.2f}", handler->GetSampleTitle(iSample), handler->GetMCHist(iSample)->Integral());

      if (handler->GetNDim(iSample) == 1) {
        handler->AddData(iSample, static_cast<TH1D*>(DUNEHists.back().get()));
      } else if (handler->GetNDim(iSample) == 2) {
        handler->AddData(iSample, static_cast<TH2D*>(DUNEHists.back().get()));
      }
    }
  }
  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

  //###############################################################################################################################
  //Lets benefit from the core code utilities 
  
  //Add samples to FitterBase
  for(auto Sample : samples){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->AddSystObj(param_handler.get());
  
  if (do_1d_llhscan) {
    MaCh3Fitter->RunLLHScan();
  }
  if (do_2d_llhscan) {
    MaCh3Fitter->Run2DLLHScan();
  }
}
