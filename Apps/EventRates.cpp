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
  auto fitMan = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create SampleHandlerFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;
  
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

  // std::vector<double> Params;
  // std::vector<double> AvgParams(480, 0.0);
  // std::vector<int> Count(480, 0);

  // TFile* Osc = TFile::Open("TrueIndChanOsc.root"); // Loading in histograms
  // TFile* Unosc = TFile::Open("TrueIndChanUnosc.root");

  // TIter next(Osc->GetListOfKeys()); // Get list of different items within oscillated data histograms (48 in total, 12 channels in 4 samples)
  // TKey* key;  // Initialise
  // int ParIndex = 0.0;

  // for(int i = 0; i < xsec->GetNumParams(); i++){ // For every param in the xsec group
  //   if(xsec->IsParFromGroup(i, "EParam")){ // If param is from our energy normalisation parameter group
  //     ParIndex = i; // Set the index of first parameter
  //     break;
  //   }
  // }

  // while ((key = (TKey*)next())) { // Go through all keys in sequence
  //   auto HistoOsc = Osc->Get<TH1D>(key->GetName()); // Getting names of histograms
  //   auto HistoUnosc = Unosc->Get<TH1D>(key->GetName()); 
  //   int NumBins = HistoOsc->GetNbinsX(); // Find number of bins (number of energy normalisation parameters for this histogram)
  //   for(int j = 1; j <= NumBins; j++) { // For each energy bin, starting from 1 to avoid the overflow bin
  //     double BinSizeOsc = HistoOsc->GetBinContent(j); // Get the number of events in specific energy bin
  //     double BinSizeUnosc = HistoUnosc->GetBinContent(j);
  //     double Param;
  //     if(BinSizeUnosc == 0) { // If no unoscillated data, set param to 0
  //       Param = 0;
  //       //xsec->SetPar(ParIndex, Param);
  //       Params.push_back(Param);
  //     }
  //     else { // If unoscillated data, calculate ratio between these as needed to induce oscillation
  //       Param = BinSizeOsc / BinSizeUnosc; 
  //       //xsec->SetPar(ParIndex, Param);
  //       Params.push_back(Param);
  //     } 
  //     //ParIndex++; // Increment parameter index to keep amending in sequence
  //   }
  // }
  // for(int p = 0; p < 1920; p++){
  //   int sample = p / 480;
  //   int channel = (p % 480) / 40;
  //   int bin = p % 40;
  //   int Index = channel * 40 + bin;
  //   AvgParams[Index] += Params[p];
  //   Count[Index] += 1;
  // }
  // for(int i = 0; i < 480; i++){
  //   AvgParams[i] /= Count[i];
  //   xsec->SetPar(ParIndex, AvgParams[i]);
  //   ParIndex++;
  // }

  // for(int k = 0; k < xsec->GetNumParams(); k++){ // For every param in the xsec group
  //   if((xsec->GetParProp(k) == 0) && xsec->IsParFromGroup(k, "EParam")){ // If it has a value of 0 and is an energy normalisation parameter
  //     xsec->ToggleFixParameter(k); // Fix these params at 0
  //   }
  // }


  std::vector<TH1*> DUNEHists;
  for(auto handler : DUNEPdfs){
    handler->Reweight();
    for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
      DUNEHists.push_back(handler->GetMCHist(iSample));

      std::string EventRateString = fmt::format("{:.2f}", handler->GetMCHist(iSample)->Integral());
      MACH3LOG_INFO("Event rate for {} : {:<5}", handler->GetSampleTitle(iSample), EventRateString);
      handler->PrintIntegral(iSample);
    }
  }

  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");
  Write1DHistogramsToFile(OutFileName, DUNEHists); 
  Write1DHistogramsToPdf(OutFileName, DUNEHists);

  //###############################################################################################################################
  //Make oscillation channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Oscillation Mode Breakdown:");
  
  //TFile* outHist = new TFile("PMNSOscDataCC.root", "recreate"); // If we want to save the individual channels from the samples, create this file

  for(auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");
      int nOscChannels = handler->GetNOscChannels(iSample);
      for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
        std::vector< KinematicCut > SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
        SelecChannel.LowerBound = iOscChan;
        SelecChannel.UpperBound = iOscChan+1;
        SelectionVec.push_back(SelecChannel);
        
        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
      
        // TString HistoName = Form("%s_%s", handler->GetSampleTitle(iSample).c_str(), handler->GetFlavourName(iSample, iOscChan).c_str());
        // Hist->SetName(HistoName); // Set the name of the histograms
        // Hist->Write(); // Save to file

        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),handler->GetFlavourName(iSample, iOscChan),Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

  //outHist->Close(); // Close our individual channel histograms

  //###############################################################################################################################
  //Make interaction channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Interaction Mode Breakdown:");

  for(auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
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

        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Modes->GetMaCh3ModeName(iModeChan),Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

  //###############################################################################################################################
}
