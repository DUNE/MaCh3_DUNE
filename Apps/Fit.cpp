#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerPDSP.h"

namespace {
void PrintCVClosureCheck(SampleHandlerFD* handler) {
  MACH3LOG_INFO("-----------------------------------------------------------------------");
  MACH3LOG_INFO("CV closure check before fitter proposals for {}", handler->GetName());
  MACH3LOG_INFO("{:<40}{:<12}{:<12}{:<12}{:<12}|",
                "Sample", "DataInt", "MCInt", "Diff", "MaxBinDiff");

  double totalData = 0.0;
  double totalMC = 0.0;
  double totalAbsDiff = 0.0;
  double globalMaxAbsDiff = 0.0;

  for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
    const auto data = handler->GetDataArray(static_cast<int>(iSample));
    const auto mc = handler->GetMCArray(static_cast<int>(iSample));

    if (data.size() != mc.size()) {
      MACH3LOG_ERROR("CV closure array size mismatch for {}: data={} MC={}",
                     handler->GetSampleTitle(iSample), data.size(), mc.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    double dataIntegral = 0.0;
    double mcIntegral = 0.0;
    double sampleAbsDiff = 0.0;
    double sampleMaxAbsDiff = 0.0;

    for (size_t iBin = 0; iBin < data.size(); ++iBin) {
      dataIntegral += data[iBin];
      mcIntegral += mc[iBin];
      const double absDiff = std::abs(mc[iBin] - data[iBin]);
      sampleAbsDiff += absDiff;
      sampleMaxAbsDiff = std::max(sampleMaxAbsDiff, absDiff);
    }

    totalData += dataIntegral;
    totalMC += mcIntegral;
    totalAbsDiff += sampleAbsDiff;
    globalMaxAbsDiff = std::max(globalMaxAbsDiff, sampleMaxAbsDiff);

    MACH3LOG_INFO("{:<40}{:<12.6f}{:<12.6f}{:<12.6f}{:<12.6f}|",
                  handler->GetSampleTitle(iSample), dataIntegral, mcIntegral,
                  mcIntegral - dataIntegral, sampleMaxAbsDiff);
  }

  MACH3LOG_INFO("{:<40}{:<12.6f}{:<12.6f}{:<12.6f}{:<12.6f}|",
                "Total", totalData, totalMC, totalMC - totalData, globalMaxAbsDiff);
  MACH3LOG_INFO("CV closure summed absolute bin difference: {:.6f}", totalAbsDiff);
  MACH3LOG_INFO("-----------------------------------------------------------------------");
}
}

int main(int argc, char * argv[]) {

  auto FitManager = MaCh3ManagerFactory(argc, argv);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  ParameterHandlerGeneric* xsec = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  for (auto handler : DUNEPdfs) {
    handler->Reweight();

    // Read real data from ROOT files (defined in SampleHandler_PDSP.yaml InputFiles)
    auto* pdsp = dynamic_cast<SampleHandlerPDSP*>(handler);
    if (pdsp == nullptr) {
      MACH3LOG_ERROR("SetupExperimentData is only implemented for SampleHandlerPDSP");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    pdsp->SetupExperimentData();
    PrintCVClosureCheck(handler);

    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
      std::string name = handler->GetSampleTitle(iSample);
      sample_names.push_back(name);
      MACH3LOG_INFO("Integrals of nominal MC hists:");
      MACH3LOG_INFO("{} : {}", name.c_str(), handler->GetMCHist(iSample)->Integral());
      MACH3LOG_INFO("--------------");
    }
  }

  if (GetFromManager(FitManager->raw()["General"]["CheckCVClosureOnly"], false)) {
    MACH3LOG_INFO("General.CheckCVClosureOnly is true; stopping before fitter construction.");
    for (auto Sample : DUNEPdfs) {
      delete Sample;
    }
    delete xsec;
    return 0;
  }
  
  //###########################################################################################################
  //MCMC

  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);
  //Start chain from random position unless continuing a chain
  if(!StartFromPreviousChain){
    if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false)) {
      xsec->ThrowParameters();
    }
  }

  //Add systematic objects
  MaCh3Fitter->AddSystObj(xsec);

  if (StartFromPreviousChain) {
    std::string PreviousChainPath = FitManager->raw()["General"]["PosFileName"].as<std::string>();
    MACH3LOG_INFO("MCMC getting starting position from: {}",PreviousChainPath);
    MaCh3Fitter->StartFromPreviousFit(PreviousChainPath);
  }
  
  //Add samples
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  
  //Run fit
  MaCh3Fitter->RunMCMC();

  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
