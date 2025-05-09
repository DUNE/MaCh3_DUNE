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

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"

void Write1DHistogramsToFile(std::string OutFileName, std::vector<TH1*> Histograms) {
  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for(auto Hist : Histograms){
    Hist->Write();
  }
  OutputFile->Close();
}

void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1*> Histograms) {
  //Remove root from end of file
  OutFileName.erase(OutFileName.find('.'));
  OutFileName+=".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName+"[").c_str());
  for(auto Hist : Histograms){
    Hist->Draw("HIST");
    c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName+"]").c_str());
}

int main(int argc, char * argv[]) {
  MaCh3Utils::MaCh3Usage(argc, argv);
  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));

  //###############################################################################################################################
  //Create samplePDFFD objects
  
  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;
  
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  //###############################################################################################################################
  //Perform reweight and print total integral

  std::vector<TH1*> DUNEHists;
  for(auto Sample : DUNEPdfs){
    Sample->reweight();
    if (Sample->GetNDim() == 1)
      DUNEHists.push_back(Sample->get1DHist());
    else if (Sample->GetNDim() == 2)
      DUNEHists.push_back(Sample->get2DHist());

    std::string EventRateString = fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetTitle(), EventRateString);

    Sample->PrintIntegral();
  }

  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");
  Write1DHistogramsToFile(OutFileName, DUNEHists); 
  Write1DHistogramsToPdf(OutFileName, DUNEHists);

  //###############################################################################################################################
  //Make oscillation channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Oscillation Mode Breakdown:");
  
  for(auto Sample : DUNEPdfs) {
    MACH3LOG_INFO("======================");
    int nOscChannels = Sample->getNMCSamples();
    for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
      std::vector< KinematicCut > SelectionVec;

      KinematicCut SelecChannel;
      SelecChannel.ParamToCutOnIt = Sample->ReturnKinematicParameterFromString("OscillationChannel");
      SelecChannel.LowerBound = iOscChan;
      SelecChannel.UpperBound = iOscChan+1;
      SelectionVec.push_back(SelecChannel);
      
      TH1* Hist = Sample->get1DVarHist(Sample->GetXBinVarName(),SelectionVec);
      MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",Sample->GetTitle(),Sample->getFlavourName(iOscChan),Hist->Integral());
    }

    TH1* Hist = Sample->get1DVarHist(Sample->GetXBinVarName());
    MACH3LOG_INFO("{:<20} : {:<20.2f}",Sample->GetTitle(),Hist->Integral());
  }

  //###############################################################################################################################
  //Make interaction channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Interaction Mode Breakdown:");

  for(auto Sample : DUNEPdfs) {
    MACH3LOG_INFO("======================");

    MaCh3Modes* Modes = Sample->GetMaCh3Modes();
    int nModeChannels = Modes->GetNModes();
    for (int iModeChan=0;iModeChan<nModeChannels;iModeChan++) {
      std::vector< KinematicCut > SelectionVec;

      KinematicCut SelecChannel;
      SelecChannel.ParamToCutOnIt = Sample->ReturnKinematicParameterFromString("Mode");
      SelecChannel.LowerBound = iModeChan;
      SelecChannel.UpperBound = iModeChan+1;
      SelectionVec.push_back(SelecChannel);

      TH1* Hist = Sample->get1DVarHist(Sample->GetXBinVarName(),SelectionVec);
      MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",Sample->GetTitle(),Modes->GetMaCh3ModeName(iModeChan),Hist->Integral());
    }

    TH1* Hist = Sample->get1DVarHist(Sample->GetXBinVarName());
    MACH3LOG_INFO("{:<20} : {:<20.2f}",Sample->GetTitle(),Hist->Integral());
  }

  //###############################################################################################################################
}
