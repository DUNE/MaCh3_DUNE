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

#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"

int main(int argc, char * argv[]) {

  auto FitManager = MaCh3ManagerFactory(argc, argv);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  //####################################################################################
  //Create sample handler + parameter_handler objects
  auto [param_handler, samples] = MaCh3DuneFactory(FitManager);

  //###########################################################################################################
  //MCMC

  auto MaCh3Fitter = std::make_unique<MR2T2>(FitManager.get());

  //Add systematic objects
  MaCh3Fitter->AddSystObj(param_handler.get());

  //Add samples
  for(auto Sample : samples){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->RunSigmaVar();

  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", M3::Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", M3::Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
