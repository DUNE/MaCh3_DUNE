#include <TH1.h>
#include "TH2.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "DUNEStyle.h"

#include <iostream>

TH1* rebinHist(TH1* hist) {
  if (std::string(hist->GetTitle()).find("BAngle_Momentum") != std::string::npos) {
    TH2D* hist2D = dynamic_cast<TH2D*>(hist);
    if (hist2D) return (TH1*)hist2D->Rebin2D(1, 2);
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

void changeAxisTitle(TAxis* axis) {
  std::string title = axis->GetTitle();
  if (title ==  "TrueQ0") axis->SetTitle("q_{0} [GeV]");
  else if (title == "TrueQ3") {
    axis->SetTitle("q_{3} [GeV/c]");
    axis->SetNdivisions(4,8,0); 
  }
  else if (title == "Particle_BAngle") axis->SetTitle("Angle to B-Field [  #circ ]");
  else if (title == "Particle_Momentum") axis->SetTitle("Momentum [GeV/c]");
  else if (title == "Particle_MomResMS") axis->SetTitle("Momentum Resolution (Multiple Scattering)");
  else if (title == "Particle_MomResYZ") axis->SetTitle("Momentum Resolution (Gluckstern)");
  else if (title == "Particle_TrackLengthYZ") axis->SetTitle("Track Length in YZ Plane (cm)");
  else if (title == "Particle_EndX") axis->SetTitle("Track End X (cm)");
  else if (title == "Particle_EndR") axis->SetTitle("Track End Radius (cm)");
}

void makeAcceptanceCorrectionPlots(const char* inputfilename) {

  TFile* inputfile = TFile::Open(inputfilename, "READ");
  if (!inputfile || inputfile->IsZombie()) {
    std::cerr << "Could not open input file " << inputfilename << std::endl;
    return;
  }

  std::string outputfilename = "Outputs/acceptance_plots/";
  std::string inputfilestring(inputfilename);
  const std::string token = "Projections";
  size_t pos = inputfilestring.find(token);
  if (pos != std::string::npos) {  // rfind==0 means "prefix at position 0"
    outputfilename += "AcceptancePlots" + inputfilestring.substr(pos + token.size());
    const std::string suffix = ".root";
    outputfilename = outputfilename.substr(0, outputfilename.size() - suffix.size()) + ".pdf";
  }
  else {
    outputfilename = "Outputs/acceptance_plots/AcceptancePlots.pdf";
  }
  const char* outputfile = strdup(outputfilename.c_str());

  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas("canvas", "Acceptance Correction Plots", 800, 600);
  gPad->SetTopMargin(0.13);
  canvas->Print(Form("%s[", outputfile));  

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
    if (obj->InheritsFrom(TH1D::Class())) rawHist->Draw("HIST");
    else rawHist->Draw("COLZ");
    rawHist->GetXaxis()->SetTitleOffset(1.3);
    canvas->Print(outputfile);
    std::cout << "Drawn histogram: " << key->GetName() << std::endl;

    if (histtitle.size()>9 && histtitle.substr(0,9) == "Accepted_") {
      int accstringpos = histname.find("Accepted_");
      std::string basename = histname.substr(0,accstringpos) + histname.substr(accstringpos+9);
      TH1* acceptedHist = (TH1*)obj;
      TH1* totalHist = (TH1*)inputfile->Get(basename.c_str());
      if (totalHist) {
        TH1* acceptedRebinned = rebinHist(acceptedHist);
        TH1* totalRebinned = rebinHist(totalHist);
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
        if (obj->InheritsFrom(TH1D::Class())) acceptanceHist->Draw("HIST");
        else acceptanceHist->Draw("COLZ");
        acceptanceHist->GetXaxis()->SetTitleOffset(1.3);
        canvas->Print(outputfile);
        std::cout << "\nDrawn acceptance histogram: " << acceptanceHist->GetName() << std::endl;
      }
    }
  }
  canvas->Print(Form("%s]", outputfile));
  inputfile->Close();
  delete inputfile;
  delete canvas;

  std::cout << "All histograms drawn in " << outputfilename << std::endl;
}
