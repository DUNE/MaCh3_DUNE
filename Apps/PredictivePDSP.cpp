// MaCh3 includes
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerPDSP.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>

int main(int argc, char *argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  const std::string PosteriorFile = Get<std::string>(
      FitManager->raw()["Predictive"]["PosteriorFile"], __FILE__, __LINE__);
  const std::string PredictiveOutputFile = GetFromManager<std::string>(
      FitManager->raw()["Predictive"]["OutputFile"],
      "PredictiveOutput.root", __FILE__, __LINE__);
  if (PredictiveOutputFile == PosteriorFile) {
    MACH3LOG_ERROR("Predictive output file '{}' must differ from posterior input file '{}'",
                   PredictiveOutputFile, PosteriorFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  const bool PriorPredictive = Get<bool>(
      FitManager->raw()["Predictive"]["PriorPredictive"], __FILE__, __LINE__);
  if (!PriorPredictive) {
    std::unique_ptr<TFile> PosteriorInput(TFile::Open(PosteriorFile.c_str(), "READ"));
    TTree* PosteriorTree = PosteriorInput && !PosteriorInput->IsZombie()
        ? PosteriorInput->Get<TTree>("posteriors")
        : nullptr;
    if (!PosteriorInput || PosteriorInput->IsZombie() ||
        PosteriorInput->TestBit(TFile::kRecovered)) {
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

  const bool UseData = GetFromManager(FitManager->raw()["General"]["Data"], false);
  const std::string AsimovTune = GetFromManager<std::string>(
      FitManager->raw()["General"]["Systematics"]["XsecAsimovTune"], "");
  if (!UseData && !AsimovTune.empty()) {
    MACH3LOG_INFO("Generating predictive Asimov data with xsec tune '{}'", AsimovTune);
    param_handler->SetTune(AsimovTune);
  }

  std::vector<TH1*> PredictionHistograms;
  for (auto handler : samples) {
    handler->Reweight();
    for (unsigned iSample = 0; iSample < handler->GetNSamples(); ++iSample) {
      const std::string name = handler->GetSampleTitle(iSample);
      TH1* DataHist = nullptr;
      if (UseData) {
        auto* PDSPHandler = dynamic_cast<SampleHandlerPDSP*>(handler);
        if (PDSPHandler == nullptr) {
          MACH3LOG_ERROR("General.Data is currently implemented for PDSP samples only");
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        DataHist = PDSPHandler->GetDataHistogramFromInputs(static_cast<int>(iSample));
      } else {
        DataHist = static_cast<TH1*>(handler->GetMCHist(iSample)->Clone((name + "_DataHist").c_str()));
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

  if (!UseData && !AsimovTune.empty()) {
    MACH3LOG_INFO("Resetting xsec parameters to PreFitValue before predictive throws");
    param_handler->SetParameters();
  }

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
