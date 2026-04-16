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

  //###############################################################################################################################
  //Create samplePDFFD objects
  auto xsec = MaCh3CovarianceFactory<ParameterHandlerGeneric>(FitManager.get(), "Xsec");
  if (CheckNodeExists(FitManager->raw(), "General", "OscillationParameters"))
  {
    auto oscpars = Get<std::vector<double>>(FitManager->raw()["General"]["OscillationParameters"], __FILE__, __LINE__);
    xsec->SetGroupOnlyParameters("Osc", oscpars);
  }

  auto DUNEPdfs = MaCh3DuneSampleFactory(FitManager, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

  MACH3LOG_INFO("=======================================================");
  for(SampleHandlerBase* handler: DUNEPdfs){
    handler->Reweight();
    for (int iSample=0;iSample<handler->GetNSamples();iSample++) {
      MACH3LOG_INFO("Event rate for {} : {:<5.2f}", handler->GetSampleTitle(iSample), handler->GetMCHist(iSample)->Integral());
    }
  }
  
  //###############################################################################################################################  
  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());
  MaCh3Fitter->AddSystObj(xsec.get());
  //Add samples
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }
  MaCh3Fitter->RunSigmaVar();
  for (auto& pdf : DUNEPdfs) {
    delete pdf;
  }

  return 0;
}
