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

  for (auto handler : DUNEPdfs) {
    handler->Reweight();

    auto* pdsp = dynamic_cast<SampleHandlerPDSP*>(handler);
    if (pdsp == nullptr) {
      MACH3LOG_ERROR("SetupExperimentData is only implemented for SampleHandlerPDSP");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    pdsp->SetupExperimentData();

    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
      const std::string name = handler->GetSampleTitle(iSample);
      MACH3LOG_INFO("Predictive data {} integral: {}", name, handler->GetDataHist(iSample)->Integral());
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
