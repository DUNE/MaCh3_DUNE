#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "tests/Comparison.h"


int main(int argc, char * argv[]) {
  MaCh3Utils::MaCh3Usage(argc, argv);

  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));
  
  bool skip_checks = false;

  if (CheckNodeExists(fitMan->raw(), "General", "Tests", "SkipChecks") ){
    skip_checks = fitMan->raw()["General"]["Tests"]["SkipChecks"].as<bool>();
  }

  //###############################################################################################################################
  //Create samplePDFFD objects
  
  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;
  
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);


  //###############################################################################################################################
  //Perform reweight and print total integral
  //###############################################################################################################################
  //Make oscillation channel breakdown
  // Initialise output file
  std::ofstream outFile("TestNewSampleOut.txt");

  outFile<<"======"<<std::endl;
  outFile<<"Total Event Rates"<<std::endl;
  outFile<<"======"<<std::endl;
  for(auto Sample : DUNEPdfs){
    Sample->reweight();
    
    std::string EventRateString = fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    outFile<<"Event rate for "<<Sample->GetName()<<": "<<EventRateString<<std::endl;;
  }


  outFile<<"======"<<std::endl;
  outFile<<"Event Rates by osc channel"<<std::endl;
  outFile<<"======"<<std::endl;

  for(auto Sample : DUNEPdfs) {
    int nOscChannels = Sample->getNMCSamples();
    for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
      std::vector< std::vector<double> > SelectionVec;

      std::vector<double> SelecChannel(3);
      SelecChannel[0] = Sample->ReturnKinematicParameterFromString("OscChannel");
      SelecChannel[1] = iOscChan;
      SelecChannel[2] = iOscChan+1;
      SelectionVec.push_back(SelecChannel);
      
      TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy",SelectionVec);
      outFile<<Sample->GetName()<<" "<<Sample->getFlavourName(iOscChan)<<": "<<Hist->Integral()<<std::endl;
    }

    TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy");
    outFile<<Sample->GetName()<<" "<<Hist->Integral()<<std::endl;
  }

  //###############################################################################################################################
  //Make interaction channel breakdown
  outFile<<"======"<<std::endl;
  outFile<<"Event Rates by interaction mode"<<std::endl;
  outFile<<"======"<<std::endl;

  for(auto Sample : DUNEPdfs) {
    int nModeChannels = kMaCh3_nModes;
    for (int iModeChan=0;iModeChan<nModeChannels;iModeChan++) {
      std::vector< std::vector<double> > SelectionVec;

      std::vector<double> SelecChannel(3);
      SelecChannel[0] = Sample->ReturnKinematicParameterFromString("Mode");
      SelecChannel[1] = iModeChan;
      SelecChannel[2] = iModeChan+1;
      SelectionVec.push_back(SelecChannel);

      TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy",SelectionVec);
      outFile<<Sample->GetName()<<" "<<MaCh3mode_ToDUNEString((MaCh3_Mode)iModeChan)<<": "<<Hist->Integral()<<std::endl;
    }

    TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy");
    outFile<<Sample->GetName()<<": "<<Hist->Integral()<<std::endl;
  }

  // Do you want to gener
  if(skip_checks){
    return 0;
  }

  std::string results_file = fitMan->raw()["General"]["Tests"]["TestResultsFile"].as<std::string>();


  // Okay now we've written we need to compare 
  bool TheSame = CompareTwoFiles(results_file, "TestNewSampleOut.txt");

  // Are we right?
  if(!TheSame) {
    MACH3LOG_CRITICAL("Different event rates detected");
    throw MaCh3Exception(__FILE__ , __LINE__ );
  } else {
    MACH3LOG_INFO("Event rates match");
  }

  return 0;
}
