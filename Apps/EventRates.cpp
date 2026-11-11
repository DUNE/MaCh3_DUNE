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

  // TFile* Osc = TFile::Open("EventRates/AllSampleTwentyTrueIndChanOsc.root"); // Loading in histograms
  // TFile* Unosc = TFile::Open("EventRates/AllSampleTwentyTrueIndChanUnosc.root");

  // TIter next(Osc->GetListOfKeys()); // Get list of different items within oscillated data histograms (48 in total, 12 channels in 4 samples)
  // TKey* key;  // Initialise

  // std::vector<std::string> KeyNames;
  // //std::vector<int> SelectedKeys = {12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
  // std::vector<int> SelectedKeys = {0, 13, 26, 39, 16, 41, 6, 31, 20, 9, 46, 35};
  //int ParIndex = 0;
  // //int KeyIndex = 0;

  // std::vector<double> NormParams(960, 0.0);
  // auto FinChain = TFile::Open("AllSampleTwentyTest.root");
  // auto Posteriors = FinChain->Get<TTree>("posteriors");
  // for (int i = 0; i < xsec->GetNumParams(); i++){
  //   std::string ChainVarName = "param_" + std::to_string(i);
  //   Posteriors->SetBranchAddress(ChainVarName.c_str(), &NormParams[i]);
  // }

  // Posteriors->GetEntry(2000234);

  // for(int i = 0; i < xsec->GetNumParams(); i++){ // For every param in the xsec group
  //   if(xsec->IsParFromGroup(i, "EParam")){ // If param is from our energy normalisation parameter group
  //     ParIndex = i; // Set the index of first parameter
  //     break;
  //   }
  // }

  // for (int i = 0; i < xsec->GetNumParams(); i++){
  //   double val = NormParams[ParIndex];
  //   xsec->SetPar(ParIndex, val);
  //   if (val == 0) xsec->ToggleFixParameter(ParIndex);
  //   ParIndex++;
  // }

  // while((key = (TKey*)next())) {
  //   KeyNames.push_back(key->GetName()); 
  // }

  // //for (int index = 0; index < 12; index++) {
  // for (int index : SelectedKeys) {
  //   const std::string& name = KeyNames[index];
  //   auto HistoOsc = Osc->Get<TH1D>(name.c_str()); // Getting names of histograms
  //   auto HistoUnosc = Unosc->Get<TH1D>(name.c_str()); 
  //   int NumBins = HistoOsc->GetNbinsX();
  //   for(int j = 1; j <= NumBins; j++) { // For each energy bin, starting from 1 to avoid the overflow bin
  //     double BinSizeOsc = HistoOsc->GetBinContent(j); // Get the number of events in specific energy bin
  //     double BinSizeUnosc = HistoUnosc->GetBinContent(j);
  //     double Param = (BinSizeUnosc == 0.0 ? 0.0 : BinSizeOsc/BinSizeUnosc);
  //     xsec->SetPar(ParIndex, Param);
  //     if (BinSizeUnosc == 0.0) xsec->ToggleFixParameter(ParIndex);
  //     ParIndex++; // Increment parameter index to keep amending in sequence
  //   }
  // }

  // for (int i = 0; i < 12; i++) {
  //   key = (TKey*) next();
  // }

  // TFile* PostFile = TFile::Open("AllSampleTwentyTest_Process.root");
  // auto HPDVal = PostFile->Get<TVectorT<double>>("Means_HPD");
  // for (int i = 0; i < HPDVal->GetNrows(); i++) {
  //   xsec->SetPar(ParIndex, (*HPDVal)[i]); 
  //   ParIndex++; 
  // }

  // while ((key = (TKey*)next())) { // Go through all keys in sequence
  //   if(KeyIndex == 24.0) break;
  //   KeyIndex++;
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
  
  TFile* outHist = new TFile("Stacked.root", "recreate"); // If we want to save the individual channels from the samples, create this file

  for(auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");
      struct HistInfo{
        TH1* hist;
        double integral;
        int oscChan;
      };
      std::vector<HistInfo> histList;
      int nOscChannels = handler->GetNOscChannels(iSample);
      for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
        std::vector< KinematicCut > SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
        SelecChannel.LowerBound = iOscChan;
        SelecChannel.UpperBound = iOscChan+1;
        SelectionVec.push_back(SelecChannel);
        
        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
        //TH1* Hist = handler->Get2DVarHist(iSample, "TrueNeutrinoEnergy", "RecoNeutrinoEnergy");
        
        TString HistoName = Form("%s_%s", handler->GetSampleTitle(iSample).c_str(), handler->GetFlavourName(iSample, iOscChan).c_str());
        Hist->SetName(HistoName); // Set the name of the histograms
        histList.push_back({Hist, Hist->Integral(), iOscChan});
        //Hist->Write(); // Save to file

        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),handler->GetFlavourName(iSample, iOscChan),Hist->Integral());
      }

      std::sort(histList.begin(), histList.end(),[](const HistInfo& a, const HistInfo& b){
        return a.integral < b.integral;
      });

      THStack* hs = new THStack(Form("hs_%d", iSample), handler->GetSampleTitle(iSample).c_str());

      for (auto& hinfo : histList){
        int iOscChan = hinfo.oscChan;
        hinfo.hist->SetLineColor(iOscChan + 2);
        if (iOscChan == 8) hinfo.hist->SetLineColor(kBlack);
        hs->Add(hinfo.hist);
      }

      TFile* Data = TFile::Open("EventRates/AllSampleTwentyOsc.root");
      TString HistName = "h" + handler->GetSampleTitle(iSample);
      TH1D* HistD = Data->Get<TH1D>(HistName);

      outHist->cd();

      TCanvas* c = new TCanvas(Form("c_%s", handler->GetSampleTitle(iSample).c_str()), "", 800, 600);
      hs->Draw("hist");
      HistD->Draw("hist same");
      HistD->SetLineStyle(2);
      HistD->SetStats(0);
      HistD->SetLineColor(kRed);
      c->BuildLegend();
      c->Write();

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

  outHist->Close(); // Close our individual channel histograms

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
