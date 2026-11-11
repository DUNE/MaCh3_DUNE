#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <cmath>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLine.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

//=====================================================================================
// Global style: call once, before making any plots. Gives you a clean, uncluttered
// look that reproduces well in print (no default ROOT grey stats box, larger fonts,
// sensible margins, thin serif-free axis labels).
//=====================================================================================
void SetThesisStyle() {
  TStyle* thesisStyle = new TStyle("Thesis", "Thesis style");

  // Canvas / pad
  thesisStyle->SetCanvasColor(kWhite);
  thesisStyle->SetPadColor(kWhite);
  thesisStyle->SetFrameBorderMode(0);
  thesisStyle->SetCanvasBorderMode(0);
  thesisStyle->SetPadBorderMode(0);
  thesisStyle->SetPadTickX(1);
  thesisStyle->SetPadTickY(1);

  // Stats / title boxes off
  thesisStyle->SetOptStat(0);
  thesisStyle->SetOptTitle(0); // we'll draw our own title with TLatex for control over placement

  // Fonts: 132 = Times-Roman, precision 2 (scales with pad, not fixed pixel size).
  // This is the standard choice for HEP theses since it visually matches the
  // Computer Modern / Times serif body text of a LaTeX document, rather than the
  // sans-serif ROOT default -- much less jarring when the plot sits in a subfigure
  // next to LaTeX-set captions.
  const int font = 132;
  thesisStyle->SetTextFont(font);
  thesisStyle->SetLabelFont(font, "XYZ");
  thesisStyle->SetTitleFont(font, "XYZ");
  thesisStyle->SetLegendFont(font);
  // Sizes bumped up a bit vs. a full-page plot: these will be scaled down when
  // placed as a subfigure, so labels need to stay legible after shrinking.
  thesisStyle->SetLabelSize(0.050, "XYZ");
  thesisStyle->SetTitleSize(0.055, "XYZ");
  thesisStyle->SetTitleOffset(1.05, "X");
  thesisStyle->SetTitleOffset(1.30, "Y");

  // Margins
  thesisStyle->SetPadTopMargin(0.08);
  thesisStyle->SetPadBottomMargin(0.14);
  thesisStyle->SetPadLeftMargin(0.14);
  thesisStyle->SetPadRightMargin(0.05);

  // Line/marker defaults
  thesisStyle->SetHistLineWidth(2);
  thesisStyle->SetEndErrorSize(0);

  thesisStyle->cd();
  gROOT->ForceStyle();
}

//=====================================================================================
// Sanitise a sample title into something safe to use as a filename: spaces and
// anything that isn't alphanumeric/underscore/hyphen become underscores.
//=====================================================================================
std::string SanitiseFileName(const std::string& in) {
  std::string out = in;
  for (auto& c : out) {
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-') c = '_';
  }
  return out;
}

//=====================================================================================
// Overlay plot for a systematic scan: every parameter value drawn on one panel.
// Sized and styled to be dropped straight into a LaTeX \subfigure -- square-ish
// canvas, Times font, no title box (caption does that job in the thesis text).
//
// Writes ONE PDF (+ one high-res PNG) PER SAMPLE, named <OutFileName>_overlay_<sample>.pdf.
//=====================================================================================
void WriteOverlayPdf(const std::string& OutFileName,
                      std::map<std::string, std::map<double, TH1*>>& SysHists,
                      const std::string& ParamName) {

  std::string BaseName = OutFileName;
  BaseName.erase(BaseName.find('.'));

  const int colours[] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1};

  for (auto& [title, valMap] : SysHists) {
    std::string PdfName = BaseName + "_overlay_" + SanitiseFileName(title) + ".pdf";
    std::string PngName = BaseName + "_overlay_" + SanitiseFileName(title) + ".png";

    // Square-ish canvas reads well as a subfigure at ~half textwidth
    auto c1 = std::unique_ptr<TCanvas>(new TCanvas(
        Form("c_overlay_%s", SanitiseFileName(title).c_str()), "c_overlay", 800, 700));
    c1->SetTopMargin(0.06);
    c1->SetBottomMargin(0.13);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->cd();

    std::vector<std::pair<double, TH1*>> sortedVals(valMap.begin(), valMap.end());
    std::sort(sortedVals.begin(), sortedVals.end(),
              [](auto& a, auto& b) { return a.first < b.first; });

    std::vector<TH1*> clones;
    double maxY = 0;
    int ic = 0;

    auto leg = std::make_unique<TLegend>(0.55, 0.66, 0.85, 0.84);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(132);
    leg->SetTextSize(0.040);

    for (auto& [val, hOrig] : sortedVals) {
      TH1* h = (TH1*)hOrig->Clone(Form("%s_clone_%.3f", hOrig->GetName(), val));
      h->SetDirectory(nullptr);

      // Divide out bin widths so differently-binned histograms are comparable
      h->Scale(1.0, "width");

      h->SetLineColor(colours[ic % 5]);
      h->SetLineWidth(3);
      h->SetMarkerColor(colours[ic % 5]);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(0.6); // small marker so the stat error bar caps are visible without cluttering the line

      // LaTeX-style axis titles: ROOT's TLatex mini-language (#it, #nu, ^{}, _{})
      // renders subscripts/italics that match what you'd typeset in the thesis body.
      h->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
      h->GetYaxis()->SetTitle("Events / Bin Width (GeV^{-1})");
      h->GetXaxis()->SetTitleSize(0.048);
      h->GetXaxis()->SetLabelSize(0.042);
      h->GetXaxis()->SetTitleOffset(1.15);
      h->GetYaxis()->SetTitleSize(0.048);
      h->GetYaxis()->SetLabelSize(0.042);
      h->GetYaxis()->SetTitleOffset(1.35);

      maxY = std::max(maxY, h->GetMaximum());
      leg->AddEntry(h, Form("%s = %.2f", ParamName.c_str(), val), "lep");
      clones.push_back(h);
      ic++;
    }

    for (size_t i = 0; i < clones.size(); ++i) {
      clones[i]->SetMaximum(maxY * 1.35);
      clones[i]->SetMinimum(0.0);
      // "HIST" draws the connecting line, "E1" overlays stat error bars with caps
      clones[i]->Draw(i == 0 ? "HIST E1" : "HIST E1 SAME");
    }
    leg->Draw();

    c1->Print(PdfName.c_str());

    // Bump resolution for the PNG export only (vector PDF is already print-quality;
    // PNG needs the extra pixels to not look soft when scaled in a document).
    c1->SetCanvasSize(1600, 1400);
    c1->Print(PngName.c_str());

    for (auto h : clones) delete h;

    MACH3LOG_INFO("Wrote overlay plot for {}: {}", title, PdfName);
  }
}

void Write1DHistogramsToFile(std::string OutFileName, std::vector<TH1*> Histograms) {
  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for (auto Hist : Histograms) {
    Hist->Write();
  }
  OutputFile->Close();
}

void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1*> Histograms) {
  OutFileName.erase(OutFileName.find('.'));
  OutFileName += ".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName + "[").c_str());
  for (auto Hist : Histograms) {
    Hist->Draw("HIST COLZ");
    c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName + "]").c_str());
}

int main(int argc, char* argv[]) {
  MaCh3Utils::MaCh3Usage(argc, argv);
  auto fitMan = MaCh3ManagerFactory(argc, argv);

  SetThesisStyle();

  //###############################################################################################################################
  //Create SampleHandlerFD objects
  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

  std::vector<TH1*> DUNEHists;
  for (auto handler : DUNEPdfs) {
    handler->Reweight();
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
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
  // MissingProtonFD systematic check: set to 0.0 and 0.2, reweight, compare
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Testing MissingProtonFD functional parameter");

  const std::string ParamName = "MissingProtonFD";
  int ParIndex = xsec->GetParIndex(ParamName); // <-- confirm this matches the real API

  if (ParIndex < 0) {
    MACH3LOG_ERROR("Could not find parameter {} in ParameterHandler, check it's defined in your systematics YAML", ParamName);
  } else {
    std::vector<double> TestValues = {0.0, 0.2};
    // sample-title -> {value -> histogram}
    std::map<std::string, std::map<double, TH1*>> SysHists;

    for (double val : TestValues) {
      MACH3LOG_INFO("------------------------------------------------------------------------");
      MACH3LOG_INFO("Setting {} = {}", ParamName, val);
      xsec->SetParProp(ParIndex, val); // <-- confirm this matches the real API

      for (auto handler : DUNEPdfs) {
        handler->Reweight();
        for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
          std::string title = handler->GetSampleTitle(iSample);
          TH1* Hist = (TH1*)handler->GetMCHist(iSample)->Clone(
              Form("%s_%s_%.2f", title.c_str(), ParamName.c_str(), val));
          Hist->SetDirectory(nullptr); // detach from any open TFile so it survives
          SysHists[title][val] = Hist;

          MACH3LOG_INFO("[{} = {:.2f}] Event rate for {:<20} : {:.2f}",
                        ParamName, val, title, Hist->Integral());
        }
      }
    }

    // Print per-bin differences so you can see the shape effect, not just the integral
    MACH3LOG_INFO("------------------------------------------------------------------------");
    MACH3LOG_INFO("Bin-by-bin comparison ({} = 0.0 vs {} = 0.2):", ParamName, ParamName);
    for (auto& [title, valMap] : SysHists) {
      TH1* h0 = valMap.at(0.0);
      TH1* h2 = valMap.at(0.2);
      MACH3LOG_INFO("Sample: {}", title);
      for (int b = 1; b <= h0->GetNbinsX(); ++b) {
        double c0 = h0->GetBinContent(b);
        double c2 = h2->GetBinContent(b);
        double diff = c2 - c0;
        MACH3LOG_INFO("  bin {:<3} : nominal(0.0)={:.3f}  shifted(0.2)={:.3f}  diff={:.3f}",
                      b, c0, c2, diff);
      }
      if (std::abs(h0->Integral() - h2->Integral()) < 1e-9 &&
          [&]() { for (int b = 1; b <= h0->GetNbinsX(); ++b) if (std::abs(h0->GetBinContent(b) - h2->GetBinContent(b)) > 1e-9) return false; return true; }()) {
        MACH3LOG_WARN("Histograms for {} are IDENTICAL between par=0.0 and par=0.2 — systematic is having NO effect!", title);
      }
    }

    // Dump both sets of histograms to file/pdf for visual inspection
    std::vector<TH1*> FlatHists;
    for (auto& [title, valMap] : SysHists)
      for (auto& [val, h] : valMap)
        FlatHists.push_back(h);

    Write1DHistogramsToFile("MissingProtonFD_test.root", FlatHists);
    Write1DHistogramsToPdf("MissingProtonFD_test.root", FlatHists);
    WriteOverlayPdf("MissingProtonFD_test.root", SysHists, ParamName);

    // Reset parameter back to nominal before continuing with the rest of the app
    xsec->SetParProp(ParIndex, 0.0);
    for (auto handler : DUNEPdfs) handler->Reweight();
  }

  //###############################################################################################################################
  //Make oscillation channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Oscillation Mode Breakdown:");
  for (auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");
      int nOscChannels = handler->GetNOscChannels(iSample);
      for (int iOscChan = 0; iOscChan < nOscChannels; iOscChan++) {
        std::vector<KinematicCut> SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("OscillationChannel");
        SelecChannel.LowerBound = iOscChan;
        SelecChannel.UpperBound = iOscChan + 1;
        SelectionVec.push_back(SelecChannel);
        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample), SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), handler->GetFlavourName(iSample, iOscChan), Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), Hist->Integral());
    }
  }

  //###############################################################################################################################
  //Make interaction channel breakdown

  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("========================================================================");
  MACH3LOG_INFO("Interaction Mode Breakdown:");

  for (auto handler : DUNEPdfs) {
    for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
      MACH3LOG_INFO("======================");

      MaCh3Modes* Modes = handler->GetMaCh3Modes();
      int nModeChannels = Modes->GetNModes();
      for (int iModeChan = 0; iModeChan < nModeChannels; iModeChan++) {
        std::vector<KinematicCut> SelectionVec;

        KinematicCut SelecChannel;
        SelecChannel.ParamToCutOnIt = handler->ReturnKinematicParameterFromString("Mode");
        SelecChannel.LowerBound = iModeChan;
        SelecChannel.UpperBound = iModeChan + 1;
        SelectionVec.push_back(SelecChannel);

        TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample), SelectionVec);
        MACH3LOG_INFO("{:<20} : {:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), Modes->GetMaCh3ModeName(iModeChan), Hist->Integral());
      }

      TH1* Hist = handler->Get1DVarHist(iSample, handler->GetXBinVarName(iSample));
      MACH3LOG_INFO("{:<20} : {:<20.2f}", handler->GetSampleTitle(iSample), Hist->Integral());
    }
  }

  //###############################################################################################################################
}