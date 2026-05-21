// MaCh3 includes
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerPDSP.h"

#include <TH1D.h>
#include <TH2D.h>

int main(int argc, char *argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  const bool UseAsimov = GetFromManager(FitManager->raw()["General"]["Asimov"], true);
  const bool UseData = GetFromManager(FitManager->raw()["General"]["Data"], false);
  if (UseAsimov == UseData) {
    MACH3LOG_ERROR("Exactly one of General.Asimov or General.Data must be true");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::vector<TH1*> PredictionHistograms;
  for (auto handler : DUNEPdfs) {
    handler->Reweight();
    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
      const std::string name = handler->GetSampleTitle(iSample);
      TH1* DataHist = nullptr;
      if (UseAsimov) {
        DataHist = static_cast<TH1*>(handler->GetMCHist(iSample)->Clone((name + "_DataHist").c_str()));
      } else {
        auto* PDSPHandler = dynamic_cast<SampleHandlerPDSP*>(handler);
        if (PDSPHandler == nullptr) {
          MACH3LOG_ERROR("General.Data is currently implemented for PDSP samples only");
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        DataHist = PDSPHandler->GetDataHistogramFromInputs(static_cast<int>(iSample));
      }
      PredictionHistograms.push_back(DataHist);

      if (handler->GetNDim(iSample) == 1) {
        handler->AddData(iSample, DataHist);
      } else if (handler->GetNDim(iSample) == 2) {
        handler->AddData(iSample, DataHist);
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
