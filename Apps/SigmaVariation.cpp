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

  //DB Sigma variations in units of each parameters Sigma
  std::vector<double> sigmaVariations = {-3, -1, 0, 1, 3};

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
  //DB Can't use the core sigma variations as it's entirely set up around the concept of multiple selections per samplePDF object
  //   Thats not the case in the FD code, which has one selection per samplePDF object
  //   Consequently have to write out own code
  
  // std::vector<std::unique_ptr<ParameterHandlerBase>> CovObjs;
  // CovObjs.emplace_back(xsec);
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
