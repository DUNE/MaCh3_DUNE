// MaCh3 includes
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"
#include "Samples/MaCh3DUNEFactory.h"

#include <TH1D.h>
#include <TH2D.h>

int main(int argc, char *argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  // For Asimov predictive plots, seed the data histograms from the
  // nominal MC prediction so the overlay has the same fake data as the fit.
  std::vector<TH1*> PredictionHistograms;
  for (auto handler : DUNEPdfs) {
    handler->Reweight();
    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
      const std::string name = handler->GetSampleTitle(iSample);
      TH1* DataHist = static_cast<TH1*>(handler->GetMCHist(iSample)->Clone((name + "_DataHist").c_str()));
      PredictionHistograms.push_back(DataHist);

      if (handler->GetNDim(iSample) == 1) {
        handler->AddData(iSample, static_cast<TH1D*>(DataHist));
      } else if (handler->GetNDim(iSample) == 2) {
        handler->AddData(iSample, static_cast<TH2D*>(DataHist));
      } else {
        MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting");
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      MACH3LOG_INFO("Predictive data seed {} integral: {}", name, DataHist->Integral());
    }
  }

  std::unique_ptr<PredictiveThrower> MaCh3Fitter = std::make_unique<PredictiveThrower>(FitManager.get());
  MaCh3Fitter->AddSystObj(xsec);
  for (auto Sample : DUNEPdfs) {
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->ProduceToys();
  MaCh3Fitter->RunPredictiveAnalysis();

  for (auto Sample : DUNEPdfs) {
    delete Sample;
  }
  delete xsec;

  return 0;
}
