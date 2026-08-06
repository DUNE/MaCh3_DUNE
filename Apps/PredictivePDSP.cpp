// MaCh3 includes
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"
#include "Samples/MaCh3DUNEFactory.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>

int main(int argc, char *argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  const std::string PosteriorFile = Get<std::string>(
      FitManager->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);
  const std::string PredictiveOutputFile = Get<std::string>(
      FitManager->raw()["General"]["OutputFile"], __FILE__, __LINE__);
  if (PredictiveOutputFile == PosteriorFile) {
    MACH3LOG_ERROR("Predictive output file '{}' must differ from posterior input file '{}'",
                   PredictiveOutputFile, PosteriorFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  const bool PriorPredictive = Get<bool>(
      FitManager->raw()["Predictive"]["PriorPredictive"], __FILE__, __LINE__);
  if (!PriorPredictive) {
    std::unique_ptr<TFile> PosteriorInput(
        M3::Open(PosteriorFile, "READ", __FILE__, __LINE__));
    TTree* PosteriorTree = PosteriorInput->Get<TTree>("posteriors");
    if (PosteriorInput->TestBit(TFile::kRecovered)) {
      MACH3LOG_ERROR("Posterior input '{}' is incomplete, still open, or required ROOT recovery",
                     PosteriorFile);
      MACH3LOG_ERROR("Wait for the fit process to finish and close the file before running predictive throws");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (PosteriorTree == nullptr || PosteriorTree->GetEntries() == 0) {
      MACH3LOG_ERROR("Posterior input '{}' does not contain a non-empty 'posteriors' tree",
                     PosteriorFile);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const unsigned int BurnInSteps = Get<unsigned int>(
        FitManager->raw()["Predictive"]["BurnInSteps"], __FILE__, __LINE__);
    const double MaximumStep = PosteriorTree->GetMaximum("step");
    if (MaximumStep <= BurnInSteps) {
      MACH3LOG_ERROR("Posterior input '{}' has maximum step {}, not beyond burn-in {}",
                     PosteriorFile, MaximumStep, BurnInSteps);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    MACH3LOG_INFO("Validated posterior input {} with {} entries",
                  PosteriorFile, PosteriorTree->GetEntries());
  }
  MACH3LOG_INFO("Writing predictive results to {}", PredictiveOutputFile);

  auto [param_handler, samples] = MaCh3DuneFactory(FitManager);

  // PredictiveThrower inherits FitterBase, which opens General.OutputFile with
  // RECREATE. Keep that output separate from the posterior chain or the chain
  // will be truncated before PredictiveThrower attempts to read it.
  std::unique_ptr<PredictiveThrower> MaCh3Fitter = std::make_unique<PredictiveThrower>(FitManager.get());
  MaCh3Fitter->AddSystObj(param_handler.get());
  for (auto Sample : samples) {
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->ProduceToys();
  MaCh3Fitter->RunPredictiveAnalysis();

  for (auto Sample : samples) {
    delete Sample;
  }

  return 0;
}
