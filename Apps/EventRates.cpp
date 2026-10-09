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

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

bool PrintIntegral = true;

void Write1DHistogramsToFile(std::string OutFileName, std::vector<std::unique_ptr<TH1>>& Histograms) {
  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for(auto& Hist : Histograms){
    Hist->Write();
  }
  OutputFile->Close();
}

void Write1DHistogramsToPdf(std::string OutFileName, std::vector<std::unique_ptr<TH1>>& Histograms) {
  //Remove root from end of file
  OutFileName.erase(OutFileName.find('.'));
  OutFileName+=".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName+"[").c_str());
  for(auto& Hist : Histograms){
    Hist->Draw("HIST");
    c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName+"]").c_str());
}

int main(int argc, char * argv[]) {
  if (argc < 2) {
    MACH3LOG_ERROR("Wrong usage of MaCh3 executable!");
    MACH3LOG_ERROR("Syntax is $: {} conf.yaml [PrintIntegral]", argv[0]);
    MACH3LOG_ERROR("Where conf.yaml is a valid config file, compatible with the manager class (manager/manager.cpp/h) and PrintIntegral is a boolean [0/1] which defines if the OscChan/Mode breakdown is printed");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  auto FitManager = std::make_unique<Manager>(M3OpenConfig(argv[1]));

  if (argc == 3) {
    PrintIntegral = std::stoi(argv[2]);
  }

  //###############################################################################################################################
  //Create SampleHandlerBase objects
  auto [param_handler, samples] = MaCh3DuneFactory(FitManager);

  //###############################################################################################################################
  //Perform reweight and print total integral

  std::vector<std::unique_ptr<TH1>> DUNEHists;
  for(auto& handler : samples){
    if (!handler){
      MACH3LOG_ERROR("Sample not set up correctly");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    handler->Reweight();
    for (int iSample=0; iSample<handler->GetNSamples(); iSample++) {
      DUNEHists.push_back(M3::Clone(handler->GetMCHist(iSample)));

      std::string EventRateString = fmt::format("{:.2f}", handler->GetMCHist(iSample)->Integral());
      MACH3LOG_INFO("Event rate for {} : {:<5}", handler->GetSampleTitle(iSample), EventRateString);
      if (PrintIntegral) {handler->PrintIntegral(iSample);}
    }
  }

  std::string OutFileName = GetFromManager<std::string>(FitManager->raw()["General"]["OutputFile"], "EventRatesOutput.root", __FILE__, __LINE__);
  Write1DHistogramsToFile(OutFileName, DUNEHists);
  Write1DHistogramsToPdf(OutFileName, DUNEHists);

  //###############################################################################################################################
  //Make oscillation channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Oscillation Mode Breakdown:");

  for(auto handler : samples) {
    for (int iSample = 0; iSample < handler->GetNSamples(); iSample++) {
      MACH3LOG_INFO("======================");
      int nOscChannels = handler->GetNOscChannels(iSample);
      for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
        std::vector< KinematicCut > SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
        SelecChannel.LowerBound = iOscChan;
        SelecChannel.UpperBound = iOscChan+1;
        SelectionVec.push_back(SelecChannel);

        auto Hist = handler->Get1DVarHist(iSample, handler->GetKinVarName(iSample, 0),SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),handler->GetFlavourName(iSample, iOscChan),Hist->Integral());
      }

      auto Hist = handler->Get1DVarHist(iSample, handler->GetKinVarName(iSample, 0));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

  //###############################################################################################################################
  //Make interaction channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Interaction Mode Breakdown:");

  for(auto handler : samples) {
    for (int iSample = 0; iSample < handler->GetNSamples(); iSample++) {
      MACH3LOG_INFO("======================");

      MaCh3Modes* Modes = handler->GetMaCh3Modes();
      int nModeChannels = Modes->GetNModes();
      for (int iModeChan=0;iModeChan<nModeChannels;iModeChan++) {
        std::vector< KinematicCut > SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("Mode");
        SelecChannel.LowerBound = iModeChan;
        SelecChannel.UpperBound = iModeChan+1;
        SelectionVec.push_back(SelecChannel);

        auto Hist = handler->Get1DVarHist(iSample, handler->GetKinVarName(iSample,0), SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Modes->GetMaCh3ModeName(iModeChan),Hist->Integral());
      }

      auto Hist = handler->Get1DVarHist(iSample, handler->GetKinVarName(iSample,0));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

  //###############################################################################################################################
}
