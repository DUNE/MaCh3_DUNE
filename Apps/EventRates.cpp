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
  //       xsec->SetPar(ParIndex, Param);
  //     }
  //     else { // If unoscillated data, calculate ratio between these as needed to induce oscillation
  //       Param = BinSizeOsc / BinSizeUnosc; 
  //       xsec->SetPar(ParIndex, Param);
  //     } 
  //     ParIndex++; // Increment parameter index to keep amending in sequence
  //   }
  // }

  // for(int k = 0; k < xsec->GetNumParams(); k++){ // For every param in the xsec group
  //   if((xsec->GetParProp(k) == 0) && xsec->IsParFromGroup(k, "EParam")){ // If it has a value of 0 and is an energy normalisation parameter
  //     xsec->ToggleFixParameter(k); // Fix these params at 0
  //   }
  // }


  std::vector<TH1*> DUNEHists;
  for(auto Sample : DUNEPdfs){
    Sample->Reweight();
    DUNEHists.push_back(Sample->GetMCHist(Sample->GetNDim()));

    std::string EventRateString = fmt::format("{:.2f}", Sample->GetMCHist(Sample->GetNDim())->Integral());
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
  
  // TFile* outHist = new TFile("TrueCCIndChanReweight.root", "recreate"); // If we want to save the individual channels from the samples, create this file

  for(auto Sample : DUNEPdfs) {
    MACH3LOG_INFO("======================");
    
    Sample->Reweight();
    int nOscChannels = Sample->GetNOscChannels();

    for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
      std::vector< KinematicCut > SelectionVec;

      KinematicCut SelecChannel;
      SelecChannel.ParamToCutOnIt = Sample->ReturnKinematicParameterFromString("OscillationChannel");
      SelecChannel.LowerBound = iOscChan;
      SelecChannel.UpperBound = iOscChan+1;
      SelectionVec.push_back(SelecChannel);
      
      TH1* Hist = Sample->Get1DVarHist(Sample->GetXBinVarName(),SelectionVec);
      TString HistoName = Form("%s_%s", Sample->GetTitle().c_str(), Sample->GetFlavourName(iOscChan).c_str());

      // Hist->SetName(HistoName); // Set the name of the histograms
      // Hist->Write(); // Save to file

      MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",Sample->GetTitle(),Sample->GetFlavourName(iOscChan),Hist->Integral());
    }

    TH1* Hist = Sample->Get1DVarHist(Sample->GetXBinVarName());
    MACH3LOG_INFO("{:<20} : {:<20.2f}",Sample->GetTitle(),Hist->Integral());
  }

  // outHist->Close(); // Close our individual channel histograms

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

      TH1* Hist = Sample->Get1DVarHist(Sample->GetXBinVarName(),SelectionVec);
      MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",Sample->GetTitle(),Modes->GetMaCh3ModeName(iModeChan),Hist->Integral());
    }

    TH1* Hist = Sample->Get1DVarHist(Sample->GetXBinVarName());
    MACH3LOG_INFO("{:<20} : {:<20.2f}",Sample->GetTitle(),Hist->Integral());
  }

  //###############################################################################################################################
}
