// MaCh3 includes
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"
#include "Samples/MaCh3DUNEFactory.h"

int main(int argc, char *argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

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
