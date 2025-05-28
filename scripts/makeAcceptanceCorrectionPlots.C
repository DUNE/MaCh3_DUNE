#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

#include <iostream>

TH1* rebinHist(TH1* hist) {
  if (std::string(hist->GetTitle()).find("BAngle_Momentum") != std::string::npos) {
    TH2D* hist2D = dynamic_cast<TH2D*>(hist);
    if (hist2D) return (TH1*)hist2D->Rebin2D(2, 4);
  }
  else if (std::string(hist->GetTitle()).find("BAngle_MomRes") != std::string::npos
    || std::string(hist->GetTitle()).find("Angular_Distribution") != std::string::npos) {
    TH2D* hist2D = dynamic_cast<TH2D*>(hist);
    if (hist2D) return (TH1*)hist2D->Rebin2D(2, 2);
  }
  else if (std::string(hist->GetTitle()).find("TrueQ3_vs_TrueQ0") != std::string::npos) {
    TH2D* hist2D = dynamic_cast<TH2D*>(hist);
    if (hist2D) return (TH1*)hist2D->Rebin2D(2, 2);
  }
  return hist;
}

void movePalette(TH1* hist) {
  gPad->Update(); // palette doesn't exist until after first draw
  TPaletteAxis* palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
  if (palette) {
    palette->SetX1NDC(0.87);  // Left edge
    palette->SetX2NDC(0.9);  // Right edge
  }
}

void changeAxisTitle(TAxis* axis) {
  std::string title = axis->GetTitle();
  if (title ==  "TrueQ0") axis->SetTitle("Q0 (GeV)");
  else if (title == "TrueQ3") axis->SetTitle("Q3 (GeV/c)");
  else if (title == "Particle_BAngle") axis->SetTitle("Angle to B-Field (  #circ )");
  else if (title == "Particle_Momentum") axis->SetTitle("Momentum (GeV)");
  else if (title == "Particle_MomResMS") axis->SetTitle("Momentum Resolution (Multiple Scattering)");
  else if (title == "Particle_MomResYZ") axis->SetTitle("Momentum Resolution (Gluckstern)");
  else if (title == "Particle_TrackLengthYZ") axis->SetTitle("Track Length in YZ Plane (cm)");
  else if (title == "Particle_EndX") axis->SetTitle("Track End X (cm)");
  else if (title == "Particle_EndR") axis->SetTitle("Track End Radius (cm)");
}

void makeAcceptanceCorrectionPlots(const char* inputfilename, const char* outputfilename = "/vols/dune/jmm224/newMaCh3/MaCh3_DUNE/outputs/AcceptancePlots.pdf") {

  TFile* inputfile = TFile::Open(inputfilename, "READ");
  if (!inputfile || inputfile->IsZombie()) {
    std::cerr << "Could not open input file " << inputfilename << std::endl;
    return;
  }

  gStyle->SetPalette(kBlueRedYellow);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.2f");
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetNumberContours(128);
  TCanvas* canvas = new TCanvas("canvas", "Acceptance Correction Plots", 900, 800);
  canvas->Print(Form("%s[", outputfilename));

  TIter next(inputfile->GetListOfKeys());
  TKey* key;

  int numKeys = inputfile->GetListOfKeys()->GetSize();
  std::cout << "Number of keys in input file: " << numKeys << std::endl;

  while ((key = (TKey*)next())) {
    std::cout << "\nFound key: " << key->GetName() << std::endl;
    TObject* obj = key->ReadObj();

    if (!obj->InheritsFrom(TH1::Class())) continue;
    TH1* rawHist = (TH1*)obj;
    std::string histtitle = rawHist->GetTitle();
    std::string histname = key->GetName();
    changeAxisTitle(rawHist->GetXaxis());
    changeAxisTitle(rawHist->GetYaxis());
    rawHist->Draw("COLZ");
    movePalette(rawHist);
    canvas->Print(outputfilename);
    std::cout << "Drawn histogram: " << key->GetName() << std::endl;

    if (histtitle.size()>9 && histtitle.substr(0,9) == "Accepted_") {
      int accstringpos = histname.find("Accepted_");
      std::string basename = histname.substr(0,accstringpos) + histname.substr(accstringpos+9);
      TH1* acceptedHist = (TH1*)obj;
      TH1* totalHist = (TH1*)inputfile->Get(basename.c_str());
      TH1* acceptedRebinned = rebinHist(acceptedHist);
      TH1* totalRebinned = rebinHist(totalHist);
      if (totalHist) {
        TH1* acceptanceHist = (TH1*)acceptedRebinned->Clone((basename+" Acceptance").c_str());
        acceptanceHist->Divide(totalRebinned);
        for(int i_x =0; i_x<acceptanceHist->GetNbinsX(); i_x++){
          for(int i_y = 0; i_y<acceptanceHist->GetNbinsY(); i_y++){
            int binnum = acceptanceHist->GetBin(i_x+1, i_y+1);
            double value = acceptanceHist->GetBinContent(binnum);
            double denominator = totalRebinned->GetBinContent(binnum);
            if (denominator == 0) acceptanceHist->SetBinContent(binnum, 0);
            else if (value == 0) acceptanceHist->SetBinContent(binnum, 1e-6); 
          }
        }
        std::string rawtitle = totalHist->GetTitle();
        acceptanceHist->SetTitle(("Acceptance_"+rawtitle).c_str());

        acceptanceHist->SetMarkerSize(0.4);
        //acceptanceHist->Draw("COLZ TEXT0");
        acceptanceHist->Draw("COLZ");
        movePalette(rawHist);
        canvas->Print(outputfilename);
        std::cout << "\nDrawn acceptance histogram: " << acceptanceHist->GetName() << std::endl;
      }
    }
  }
  canvas->Print(Form("%s]", outputfilename));
  inputfile->Close();
  delete inputfile;
  delete canvas;

  std::cout << "All histograms drawn in " << outputfilename << std::endl;
}
