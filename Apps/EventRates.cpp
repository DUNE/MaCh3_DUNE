// #include <iostream>
// #include <chrono>
// #include <iomanip>
// #include <vector>

// #include <TH1D.h>
// #include <THStack.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TRint.h>
// #include <TLegend.h>
// #include <TColor.h>
// #include <TMath.h>

// #include "Samples/MaCh3DUNEFactory.h"
// #include "Samples/StructsDUNE.h"
// #include "Fitters/MaCh3Factory.h"

// void WriteOverlayPdf(const std::string& OutFileName,
//                       std::map<std::string, std::map<double, TH1*>>& SysHists,
//                       const std::string& ParamName) {

//   std::string PdfName = OutFileName;
//   PdfName.erase(PdfName.find('.'));
//   PdfName += "_overlay.pdf";

//   gStyle->SetOptStat(0);

//   auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c_overlay", "c_overlay", 900, 700));
//   c1->cd();
//   c1->Print(std::string(PdfName + "[").c_str());

//   for (auto& [title, valMap] : SysHists) {
//     // Clone so we don't mutate the histograms other code might still use
//     std::vector<std::pair<double, TH1*>> sortedVals(valMap.begin(), valMap.end());
//     std::sort(sortedVals.begin(), sortedVals.end(),
//               [](auto& a, auto& b){ return a.first < b.first; });

//     std::vector<TH1*> clones;
//     double maxY = 0;

//     int colours[] = {kBlack, kRed+1, kAzure+2, kGreen+2, kMagenta+1};
//     int ic = 0;

//     auto leg = std::make_unique<TLegend>(0.65, 0.7, 0.88, 0.88);
//     leg->SetBorderSize(0);
//     leg->SetFillStyle(0);

//     for (auto& [val, hOrig] : sortedVals) {
//       TH1* h = (TH1*)hOrig->Clone(Form("%s_clone_%.3f", hOrig->GetName(), val));
//       h->SetDirectory(nullptr);

//       // Divide out bin widths: content -> content / bin_width
//       h->Scale(1.0, "width");

//       h->SetLineColor(colours[ic % 5]);
//       h->SetLineWidth(2);
//       h->SetMarkerColor(colours[ic % 5]);
//       h->SetMarkerStyle(20);
//       h->SetTitle(Form("%s;%s;Events / bin width", title.c_str(), h->GetXaxis()->GetTitle()));

//       maxY = std::max(maxY, h->GetMaximum());
//       leg->AddEntry(h, Form("%s = %.2f", ParamName.c_str(), val), "l");
//       clones.push_back(h);
//       ic++;
//     }

//     for (size_t i = 0; i < clones.size(); ++i) {
//       clones[i]->SetMaximum(maxY * 1.3);
//       clones[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
//     }
//     leg->Draw();
//     c1->Print(PdfName.c_str());

//     for (auto h : clones) delete h;
//   }

//   c1->Print(std::string(PdfName + "]").c_str());
//   MACH3LOG_INFO("Wrote overlay PDF: {}", PdfName);
// }

// void Write1DHistogramsToFile(std::string OutFileName, std::vector<TH1*> Histograms) {
//   auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutFileName.c_str(), "RECREATE"));
//   OutputFile->cd();
//   for(auto Hist : Histograms){
//     Hist->Write();
//   }
//   OutputFile->Close();
// }

// void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1*> Histograms) {
//   //Remove root from end of file
//   OutFileName.erase(OutFileName.find('.'));
//   OutFileName+=".pdf";

//   auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
//   c1->cd();
//   c1->Print(std::string(OutFileName+"[").c_str());
//   for(auto Hist : Histograms){
//     Hist->Draw("HIST COLZ");
//     c1->Print(OutFileName.c_str());
//   }
//   c1->Print(std::string(OutFileName+"]").c_str());
// }

// int main(int argc, char * argv[]) {
//   MaCh3Utils::MaCh3Usage(argc, argv);
//   auto fitMan = MaCh3ManagerFactory(argc, argv);

//   //###############################################################################################################################
//   //Create SampleHandlerFD objects
  
//   ParameterHandlerGeneric* xsec = nullptr;
  
//   std::vector<SampleHandlerFD*> DUNEPdfs;
//   MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

//   //###############################################################################################################################
//   //Perform reweight and print total integral

//   std::vector<TH1*> DUNEHists;
//   for(auto handler : DUNEPdfs){
//     handler->Reweight();
//     for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
//       DUNEHists.push_back(handler->GetMCHist(iSample));

//       std::string EventRateString = fmt::format("{:.2f}", handler->GetMCHist(iSample)->Integral());
//       MACH3LOG_INFO("Event rate for {} : {:<5}", handler->GetSampleTitle(iSample), EventRateString);
//       handler->PrintIntegral(iSample);
//     }
//   }

//   std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");
//   Write1DHistogramsToFile(OutFileName, DUNEHists); 
//   Write1DHistogramsToPdf(OutFileName, DUNEHists);


//    //###############################################################################################################################
//   // MissingProtonFD systematic check: set to 0.0 and 0.2, reweight, compare
//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("Testing MissingProtonFD functional parameter");

//   const std::string ParamName = "MissingProtonFD";
//   int ParIndex = xsec->GetParIndex(ParamName); // <-- confirm this matches the real API

//   if (ParIndex < 0) {
//     MACH3LOG_ERROR("Could not find parameter {} in ParameterHandler, check it's defined in your systematics YAML", ParamName);
//   } else {
//     std::vector<double> TestValues = {0.0, 0.2};
//     // sample-title -> {value -> histogram}
//     std::map<std::string, std::map<double, TH1*>> SysHists;

//     for (double val : TestValues) {
//       MACH3LOG_INFO("------------------------------------------------------------------------");
//       MACH3LOG_INFO("Setting {} = {}", ParamName, val);
//       xsec->SetParProp(ParIndex, val); // <-- confirm this matches the real API

//       for (auto handler : DUNEPdfs) {
//         handler->Reweight();
//         for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
//           std::string title = handler->GetSampleTitle(iSample);
//           TH1* Hist = (TH1*)handler->GetMCHist(iSample)->Clone(
//               Form("%s_%s_%.2f", title.c_str(), ParamName.c_str(), val));
//           Hist->SetDirectory(nullptr); // detach from any open TFile so it survives
//           SysHists[title][val] = Hist;

//           MACH3LOG_INFO("[{} = {:.2f}] Event rate for {:<20} : {:.2f}",
//                         ParamName, val, title, Hist->Integral());
//         }
//       }
//     }

//     // Print per-bin differences so you can see the shape effect, not just the integral
//     MACH3LOG_INFO("------------------------------------------------------------------------");
//     MACH3LOG_INFO("Bin-by-bin comparison ({} = 0.0 vs {} = 0.2):", ParamName, ParamName);
//     for (auto& [title, valMap] : SysHists) {
//       TH1* h0 = valMap.at(0.0);
//       TH1* h2 = valMap.at(0.2);
//       MACH3LOG_INFO("Sample: {}", title);
//       for (int b = 1; b <= h0->GetNbinsX(); ++b) {
//         double c0 = h0->GetBinContent(b);
//         double c2 = h2->GetBinContent(b);
//         double diff = c2 - c0;
//         MACH3LOG_INFO("  bin {:<3} : nominal(0.0)={:.3f}  shifted(0.2)={:.3f}  diff={:.3f}",
//                       b, c0, c2, diff);
//       }
//       if (std::abs(h0->Integral() - h2->Integral()) < 1e-9 &&
//           [&](){ for(int b=1;b<=h0->GetNbinsX();++b) if(std::abs(h0->GetBinContent(b)-h2->GetBinContent(b))>1e-9) return false; return true; }()) {
//         MACH3LOG_WARN("Histograms for {} are IDENTICAL between par=0.0 and par=0.2 — systematic is having NO effect!", title);
//       }
//     }

//     // Dump both sets of histograms to file/pdf for visual inspection
//     std::vector<TH1*> FlatHists;
//     for (auto& [title, valMap] : SysHists)
//       for (auto& [val, h] : valMap)
//         FlatHists.push_back(h);

//     Write1DHistogramsToFile("MissingProtonFD_test.root", FlatHists);
//     Write1DHistogramsToPdf("MissingProtonFD_test.root", FlatHists);
//     WriteOverlayPdf("MissingProtonFD_test.root", SysHists, ParamName);
//     // Reset parameter back to nominal before continuing with the rest of the app
//     xsec->SetParProp(ParIndex, 0.0);
//     for (auto handler : DUNEPdfs) handler->Reweight();
//   }

//   //###############################################################################################################################
//   //Make oscillation channel breakdown

//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("Oscillation Mode Breakdown:");
  
//   for(auto handler : DUNEPdfs) {
//     for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
//       MACH3LOG_INFO("======================");
//       int nOscChannels = handler->GetNOscChannels(iSample);
//       for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
//         std::vector< KinematicCut > SelectionVec;

//         KinematicCut SelecChannel;
//         SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
//         SelecChannel.LowerBound = iOscChan;
//         SelecChannel.UpperBound = iOscChan+1;
//         SelectionVec.push_back(SelecChannel);
        
//         TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
//         MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),handler->GetFlavourName(iSample, iOscChan),Hist->Integral());
//       }

//       TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
//       MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
//     }
//   }

//   //###############################################################################################################################
//   //Make interaction channel breakdown

//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("========================================================================");
//   MACH3LOG_INFO("Interaction Mode Breakdown:");

//   for(auto handler : DUNEPdfs) {
//     for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
//       MACH3LOG_INFO("======================");

//       MaCh3Modes* Modes = handler->GetMaCh3Modes();
//       int nModeChannels = Modes->GetNModes();
//       for (int iModeChan=0;iModeChan<nModeChannels;iModeChan++) {
//         std::vector< KinematicCut > SelectionVec;

//         KinematicCut SelecChannel;
//         SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("Mode");
//         SelecChannel.LowerBound = iModeChan;
//         SelecChannel.UpperBound = iModeChan+1;
//         SelectionVec.push_back(SelecChannel);

//         TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample),SelectionVec);
//         MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Modes->GetMaCh3ModeName(iModeChan),Hist->Integral());
//       }

//       TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
//       MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
//     }
//   }

//   //###############################################################################################################################
// }
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
#include "Samples/SampleHandlerBeamOffAxis.h"   

void Write1DHistogramsToFile(std::string OutFileName,
                             std::vector<TH1*> Histograms) {

  std::cout << "Writing histograms to: " << OutFileName << std::endl;

  TFile* OutputFile = TFile::Open(OutFileName.c_str(), "RECREATE");

  std::cout << "TFile pointer: " << OutputFile << std::endl;

  if (!OutputFile) {
    std::cerr << "ERROR: TFile::Open returned nullptr" << std::endl;
    return;
  }

  if (OutputFile->IsZombie()) {
    std::cerr << "ERROR: Output file is a zombie" << std::endl;
    delete OutputFile;
    return;
  }

  std::cout << "File opened successfully" << std::endl;

  OutputFile->cd();

  std::cout << "cd() successful" << std::endl;

  for (auto Hist : Histograms) {
    if (!Hist) {
      std::cerr << "ERROR: null histogram!" << std::endl;
      continue;
    }

    std::cout << "Writing " << Hist->GetName() << std::endl;
    Hist->Write();
  }

  OutputFile->Close();
  delete OutputFile;
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

  std::vector<TH1*> DUNEHists;
  for(auto handler : DUNEPdfs){
    handler->Reweight();
    for (int iSample=0; iSample<handler->GetNsamples(); iSample++) {
      DUNEHists.push_back(handler->GetMCHist(iSample));

      std::string EventRateString = fmt::format("{:.2f}", handler->GetMCHist(iSample)->Integral());
      MACH3LOG_INFO("Event rate for {} : {:<5}", handler->GetSampleTitle(iSample), EventRateString);
      handler->PrintIntegral(iSample);
    }
  //   if (auto* offAxisHandler = dynamic_cast<dune::beamoffaxis::SampleHandlerBeamOffAxis*>(handler)) {
  //   offAxisHandler->MakeFancyPlot();
  // }
  }

  

  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");
  Write1DHistogramsToFile(OutFileName, DUNEHists); 
  Write1DHistogramsToPdf(OutFileName, DUNEHists);

  //###############################################################################################################################
  //Make oscillation channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Oscillation Mode Breakdown:");
  
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
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),handler->GetFlavourName(iSample, iOscChan),Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}",handler->GetSampleTitle(iSample),Hist->Integral());
    }
  }

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