#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TClass.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

#include <iostream>
#include <string>
#include <map>
#include <regex>
#include <utility>

int badfracmetric(double acceptance_threshold = 0.1) {
  // Open input file - should just contain histograms which only differ by whether they are accepted, and their enu range
  TFile* inputfile = TFile::Open("Projections_CC_BField0_5_FidRad_160_FidLen_209_InstRad_249_45_InstLen_259.root", "READ");
  if (inputfile->IsZombie()) { 
    std::cerr << "Could not open file!\n"; 
    return 1; 
  }

  // Form of enu range to search for in histogram titles
  std::regex enuRegex(R"(\(([-+]?\d*\.?\d+)\s*<=\s*Enu\s*<\s*([-+]?\d*\.?\d+)\))");

  // Map from energy range (x,y) → pair of histogram names (raw, accepted)
  std::map<std::pair<double,double>, std::pair<std::string,std::string>> enuMap;

  // Loop through hists in file to create the map
  TIter nextKey(inputfile->GetListOfKeys());
  while (TKey* key = (TKey*)nextKey()) {
    // Get histogram
    TObject* obj = key->ReadObj();
    if (!obj->InheritsFrom(TH2D::Class())) { delete obj; continue; }
    TH2D* hist = static_cast<TH2D*>(obj);
    std::string title = hist->GetTitle();
    std::string name  = hist->GetName();

    // Check for an enu range in the title
    std::smatch match;
    if (!std::regex_search(title, match, enuRegex)) {
      std::cerr << "Skipping histogram with no Enu range: " << title << std::endl;
      delete obj; 
      continue;
    }

    double enu_lower = std::stod(match[1]);
    double enu_upper = std::stod(match[2]);
    auto keyRange = std::make_pair(enu_lower, enu_upper);

    auto it = enuMap.find(keyRange);
    if (it == enuMap.end()) {
      enuMap[keyRange] = std::make_pair(name, "");
    } else {
      if (!it->second.second.empty()) {
        std::cerr << "Error: more than two histograms for range (" << enu_lower << ", " << enu_upper << ")" << std::endl;
        return 1;
      }
      it->second.second = name;
      if (name.find("Accepted") == std::string::npos) {
        if (it->second.first.find("Accepted")) {
          std::swap(it->second.first, it->second.second);
        } 
        else {
          std::cerr << "Neither hist in range (" << enu_lower << ", " << enu_upper << ") is accepted." << std::endl;
          return 1;
        }
      } 
    }
    delete obj;
  }

  // TGraphAsymmErrors parameters
  std::vector<double> x, y, exl, exh, eyl, eyh;

  // Loop through map
  for (auto& [range, names] : enuMap) {
    auto [enu_low, enu_high] = range;
    auto [rawhistname, acchistname] = names;

    // Get hists
    TH2D* rawhist = static_cast<TH2D*>(inputfile->Get(rawhistname.c_str()));
    TH2D* acchist = static_cast<TH2D*>(inputfile->Get(acchistname.c_str()));
    int nxbins, nybins;
    if (rawhist->GetNbinsX() == acchist->GetNbinsX() && rawhist->GetNbinsY() == acchist->GetNbinsY()) {
      nxbins = rawhist->GetNbinsX();
      nybins = rawhist->GetNbinsY();
    }
    else {
      std::cerr << rawhistname << " and " << acchistname << " have different binnings" << std::endl;
      return 1;
    }

    double num_events = 0;
    double num_bad_events = 0;

    for (int xbin=1; xbin<=nxbins; xbin++) {
      for (int ybin=1; ybin<=nybins; ybin++) {
        double n_raw = rawhist->GetBinContent(xbin, ybin);
        double n_acc = acchist->GetBinContent(xbin, ybin);
        if (n_raw == 0) continue; 
        num_events += n_raw;

        double acceptance = n_acc / n_raw;
        if (acceptance < acceptance_threshold) num_bad_events += n_raw;
      }
    }

    x.push_back((enu_high + enu_low) / 2);
    exl.push_back((enu_high - enu_low) / 2);
    exh.push_back((enu_high - enu_low) / 2);
    y.push_back(num_bad_events / num_events);
    eyl.push_back(0.);
    eyh.push_back(0.);
  }

  // Create output file and canvas
  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas("canvas", "Acceptance Correction Plots", 800, 600);
  gPad->SetTopMargin(0.13);
  const char* outputfilename = "badfracmetric.pdf";
  canvas->Print(Form("%s[", outputfilename));

  // Draw graph points
  TGraphAsymmErrors* fracgraph = new TGraphAsymmErrors(x.size(), x.data(), y.data(), exl.data(), exh.data(), eyl.data(), eyh.data());
  fracgraph->SetLineColor(kCyan-2);
  fracgraph->SetLineWidth(1);
  fracgraph->Draw();

  // Draw Phase II stat uncertainty line
  const int npoints = 3;
  double xpoints[npoints] = {0.0, 2.5, 5.0};
  double ypoints[npoints] = {0.03, 0.03, 0.08}; // flat then linear rise
  
  auto stat_legend = new TLegend(0.2, 0.7, 0.54, 0.82);
  TGraph* line = new TGraph(npoints, xpoints, ypoints);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  stat_legend->AddEntry(line, "Phase II stat uncertainty", "l");
  line->Draw("L SAME");
  stat_legend->Draw();

  canvas->Print(outputfilename);

  // Close everything
  canvas->Print(Form("%s]", outputfilename));
  delete canvas;
  inputfile->Close();
  delete inputfile;

  return 0;
}

