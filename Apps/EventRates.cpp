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

  TFile* Osc = TFile::Open("TrueIndChanOsc.root");
  TFile* Unosc = TFile::Open("TrueIndChanUnosc.root");

  TIter next(Osc->GetListOfKeys()); 
  TKey* key; 
  int ParIndex = 0.0; // 293.0 if all params included
  //int KeyCount = 0.0;

  for(int i = 0; i < xsec->GetNumParams(); i++){
    if( xsec->IsParFromGroup(i, "EParam")){
      ParIndex = i;
      break;
    }
  }

  while ((key = (TKey*)next())) { 
    //if (KeyCount >= 12) break;

    auto HistoOsc = Osc->Get<TH1D>(key->GetName());
    auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

    int NumBins = HistoOsc->GetNbinsX();

    for(int j = 1; j <= NumBins; j++) {
      double BinSizeOsc = HistoOsc->GetBinContent(j);
      double BinSizeUnosc = HistoUnosc->GetBinContent(j);

      double Param;
      if(BinSizeUnosc == 0) {
        Param = 0;
        xsec->SetPar(ParIndex, Param);
      }
      else {
        Param = BinSizeOsc / BinSizeUnosc; 
        xsec->SetPar(ParIndex, Param);
      }
        
      ParIndex++;
    }

    //KeyCount++;

  }


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
  
  // Include for producing plots with variable param values
  
  // TFile* Osc = TFile::Open("IndChanOsc.root");
  // TFile* Unosc = TFile::Open("IndChanUnosc.root");

  // TIter next(Osc->GetListOfKeys()); 
  // TKey* key; 
  // int ParIndex = 0.0; // 293.0 if all params included
  // int KeyCount = 0.0;

  // while ((key = (TKey*)next())) { 
  //   if (KeyCount >= 12) break;

  //   auto HistoOsc = Osc->Get<TH1D>(key->GetName());
  //   auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

  //   int NumBins = HistoOsc->GetNbinsX();

  //   for(int j = 1; j <= NumBins; j++) {
  //     double BinSizeOsc = HistoOsc->GetBinContent(j);
  //     double BinSizeUnosc = HistoUnosc->GetBinContent(j);

  //     double Param;
  //     if(BinSizeUnosc == 0) {
  //       Param = 0;
  //       xsec->SetPar(ParIndex, Param);
  //     }
  //     else {
  //       Param = BinSizeOsc / BinSizeUnosc; 
  //       xsec->SetPar(ParIndex, Param);
  //     }
        
  //     ParIndex++;
  //   }

  //   KeyCount++;

  // }
  
  // TFile* outHist = new TFile("TrueCCIndChanReweight.root", "recreate"); 

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

      // Hist->SetName(HistoName);
      // Hist->Write();

      MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",Sample->GetTitle(),Sample->GetFlavourName(iOscChan),Hist->Integral());
    }

    TH1* Hist = Sample->Get1DVarHist(Sample->GetXBinVarName());
    MACH3LOG_INFO("{:<20} : {:<20.2f}",Sample->GetTitle(),Hist->Integral());
  }

  // outHist->Close();

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
