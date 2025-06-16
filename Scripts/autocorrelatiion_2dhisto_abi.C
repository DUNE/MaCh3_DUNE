#include "TFile.h"
#include "TDirectoryFile.h"
#include "TKey.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include <iostream>
#include <vector>
#include <set>
#include <yaml-cpp/yaml.h> 

struct BinDef {
  int index;
  double q0_min, q0_max;
  double q3_min, q3_max;
};

void extract_q0q3_bins_from_yaml(const std::string& yaml_file,
  std::vector<BinDef>& binDefs,
  std::vector<double>& q0_edges,
  std::vector<double>& q3_edges) {
  std::set<double> q0_set, q3_set;
  binDefs.clear();

  try {
    YAML::Node config = YAML::LoadFile(yaml_file);
    const auto& systs = config["Systematics"];
    int index = 0;

    for (const auto& systematic : systs) {
      const auto& sys = systematic["Systematic"];
      const auto& cuts = sys["KinematicCuts"];

      double q0_min = -1, q0_max = -1, q3_min = -1, q3_max = -1;

      for (const auto& cut : cuts) {
        if (cut["q0"]) {
          const auto& q0 = cut["q0"];
          if (q0.IsSequence() && q0.size() == 2) {
            q0_min = q0[0].as<double>();
            q0_max = q0[1].as<double>();
            q0_set.insert(q0_min);
            q0_set.insert(q0_max);
          }
        }
        if (cut["q3"]) {
          const auto& q3 = cut["q3"];
          if (q3.IsSequence() && q3.size() == 2) {
            q3_min = q3[0].as<double>();
            q3_max = q3[1].as<double>();
            q3_set.insert(q3_min);
            q3_set.insert(q3_max);
          }
        }
      }

      if (q0_min >= 0 && q3_min >= 0) {
        binDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
      }
    }

    q0_edges.assign(q0_set.begin(), q0_set.end());
    q3_edges.assign(q3_set.begin(), q3_set.end());

    std::cout << " Parsed " << binDefs.size() << " bins from YAML\n";
    std::cout << "   q0 bins: " << q0_edges.size() - 1
              << ", q3 bins: " << q3_edges.size() - 1 << "\n";
  } catch (const std::exception& e) {
    std::cerr << " YAML parsing failed: " << e.what() << "\n";
  }
}

void autocorrelatiion_2dhisto_abi(const TString& diagfile, const TString& yaml_file, const TString& output) {
  std::vector<BinDef> binDefs;
  std::vector<double> q0_edges, q3_edges;

  extract_q0q3_bins_from_yaml(yaml_file.Data(), binDefs, q0_edges, q3_edges);

  int nbins_q0 = q0_edges.size() - 1;
  int nbins_q3 = q3_edges.size() - 1;

  TH2D* h_integrals = new TH2D("h_integrals", "Autocorrelation Integrals; q_{3} [GeV]; q_{0} [GeV]",
                               nbins_q3, &q3_edges[0], nbins_q0, &q0_edges[0]);

  TFile *fin = new TFile(diagfile, "OPEN");
  if (fin->IsZombie()) {
    std::cerr << "Error opening file " << diagfile << std::endl;
    return;
  }

  TDirectoryFile* autocor = (TDirectoryFile*)fin->Get("Auto_corr");
  if (!autocor) {
    std::cerr << "Could not find directory Auto_corr" << std::endl;
    return;
  }

  for (const BinDef& bin : binDefs) {
    TString histname = Form("xsec_%d", bin.index);
    TH1D* hist = dynamic_cast<TH1D*>(autocor->Get(histname));
    if (!hist) {
      std::cerr << "⚠️ Skipping missing hist: " << histname << std::endl;
      continue;
    }

    double integral = hist->Integral();

    // Find bin centers
    double q0_center = 0.5 * (bin.q0_min + bin.q0_max);
    double q3_center = 0.5 * (bin.q3_min + bin.q3_max);

    h_integrals->Fill(q3_center, q0_center, integral);
  }

  TCanvas* c = new TCanvas("c", "Autocorrelation Integrals", 800, 600);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  h_integrals->Draw("COLZ TEXT");
  c->SaveAs(output + ".png");

  fin->Close();
}
