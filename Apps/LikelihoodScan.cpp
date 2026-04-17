#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <TH3D.h>
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
  //Create samplePDFFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;
  
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight, print total integral, and set data

  std::vector<TH1*> DUNEHists;
  for(auto handler : DUNEPdfs){
    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
      handler->Reweight();
      DUNEHists.push_back(handler->GetMCHist(iSample));
      MACH3LOG_INFO("Event rate for {} : {:<5.2f}", handler->GetSampleTitle(iSample), handler->GetMCHist(iSample)->Integral());

      if (handler->GetNDim(iSample) == 1) {
        handler->AddData(iSample, (TH1D*)DUNEHists.back());
      } else if (handler->GetNDim(iSample) == 2) {
        handler->AddData(iSample, (TH2D*)DUNEHists.back());
      } else if (handler->GetNDim(iSample) == 3) {
        handler->AddData(iSample, (TH3D*)DUNEHists.back());
      }
    }
  }
  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

  //###############################################################################################################################
  //Lets benefit from the core code utilities 
  
  //Add samples to FitterBase
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->AddSystObj(xsec);
  
  if (do_1d_llhscan) {
    MaCh3Fitter->RunLLHScan();
  }
  if (do_2d_llhscan) {
    MaCh3Fitter->Run2DLLHScan();
  }
}
