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

#include "Fitters/mcmc.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"

int main(int argc, char * argv[]) {
  if(argc == 1){
    MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
    return 1;
  }
  manager* FitManager = new manager(argv[1]);

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
  ParameterHandlerOsc* osc = nullptr;
  
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);

  //###############################################################################################################################
  //Perform reweight, print total integral, and set data
  std::vector<double> oscpars = FitManager->raw()["General"]["OscillationParameters"].as<std::vector<double>>();
  for (int i = 0; i < oscpars.size(); i++)
    osc->SetPar(i, oscpars.at(i));

  std::vector<TH1*> DUNEHists;
  for(auto Sample : DUNEPdfs){
    Sample->Reweight();
    if (Sample->GetNDim() == 1) {
      DUNEHists.push_back(Sample->Get1DHist());
    } else if (Sample->GetNDim() == 2) {
      DUNEHists.push_back(Sample->Get2DHist());
    }

    MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetTitle(), Sample->Get1DHist()->Integral());
    if (Sample->GetNDim() == 1) {
      Sample->AddData((TH1D*)DUNEHists.back());
    } else if (Sample->GetNDim() == 2) {
      Sample->AddData((TH2D*)DUNEHists.back());
    }
  }
  std::unique_ptr<FitterBase> MaCh3Fitter = std::make_unique<mcmc>(FitManager);

  //###############################################################################################################################
  //Lets benefit from the core code utilities 
  
  //Add samples to FitterBase
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->AddSystObj(osc);
  MaCh3Fitter->AddSystObj(xsec);
  
  if (do_1d_llhscan) {
    MaCh3Fitter->RunLLHScan();
  }
  if (do_2d_llhscan) {
    MaCh3Fitter->Run2DLLHScan();
  }
}
