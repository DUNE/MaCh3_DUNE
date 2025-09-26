// MaCh3 includes
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"
#include "Samples/MaCh3DUNEFactory.h"

int main(int argc, char *argv[]) {
  // Initialise manger responsible for config handling
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  ParameterHandlerGeneric* xsec = nullptr;

  std::vector<SampleHandlerFD*> mySamples;
  MakeMaCh3DuneInstance(FitManager.get(), mySamples, xsec);
  // Create MCMC Class
  std::unique_ptr<PredictiveThrower> MaCh3Fitter = std::make_unique<PredictiveThrower>(FitManager.get());
  // Add covariance to PredictiveThrower
  MaCh3Fitter->AddSystObj(xsec);
  for (auto Sample: mySamples) {
    MaCh3Fitter->AddSampleHandler(Sample);
  }
  // Run Predictive analysis
  MaCh3Fitter->ProduceToys();
  MaCh3Fitter->RunPredictiveAnalysis();

  for (auto Sample: mySamples) {
    delete Sample;
  }
  return 0;
}