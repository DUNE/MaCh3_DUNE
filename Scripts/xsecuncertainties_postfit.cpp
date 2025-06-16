#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TROOT.h>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " chain_output.root [optional h_q0q3.root]" << std::endl;
    return 1;
  }

  TFile* chainFile = new TFile(argv[1]);
  TTree* post = (TTree*)chainFile->Get("posteriors");
  if (!post) {
    std::cerr << "ERROR: Could not find 'posteriors' TTree in file " << argv[1] << std::endl;
    return 1;
  }

  // === Set binning (adjust to your analysis) ===
  std::vector<double> q0_edges, q3_edges;
  for (int i = 0; i <= 50; ++i) {
    q0_edges.push_back(0.1 * i);  // 0 to 5.0
    q3_edges.push_back(0.1 * i);
  }

  int nbins_q0 = q0_edges.size() - 1;
  int nbins_q3 = q3_edges.size() - 1;
  int nParams = 1275; //nbins_q0 * nbins_q3;

  // === Initialize histograms ===
  TH2D* h_param_mean = new TH2D("h_param_mean", "Posterior Mean; q_{3} [GeV]; q_{0} [GeV]",
                                nbins_q3, &q3_edges[0], nbins_q0, &q0_edges[0]);
  TH2D* h_param_stddev = (TH2D*)h_param_mean->Clone("h_param_stddev");
  h_param_stddev->SetTitle("Posterior Std. Dev.");

  // Get binning
int nbins_q0 = h_q0q3->GetNbinsY();
int nbins_q3 = h_q0q3->GetNbinsX();
int nParams = xsec->GetNumParams();

// === Histogram: 1. Event Rate (already created as h_q0q3) ===
TH2D* h_event_rate = (TH2D*)h_q0q3->Clone("h_event_rate");
h_event_rate->SetTitle("1. Nominal Event Rate; q_{3} [GeV]; q_{0} [GeV]");

// === Histogram: 2. Initial Parameter Values (before freezing) ===
TH2D* h_param_init = (TH2D*)h_q0q3->Clone("h_param_init");
h_param_init->Reset();
h_param_init->SetTitle("2. Parameter Values Before Freezing");

for (int i = 0; i < nParams; ++i) {
  int bin_q0 = i / nbins_q3 + 1;
  int bin_q3 = i % nbins_q3 + 1;
  if (bin_q0 > nbins_q0) continue;

  double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
  double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
  if (q0 > q3) continue;

  h_param_init->SetBinContent(bin_q3, bin_q0, xsec->getParameterValue(i));
}

// === Histogram: 3. Parameter Values After Freezing ===
TH2D* h_param_frozen = (TH2D*)h_q0q3->Clone("h_param_frozen");
h_param_frozen->Reset();
h_param_frozen->SetTitle("3. Parameter Values After Freezing");

for (int i = 0; i < nParams; ++i) {
  int bin_q0 = i / nbins_q3 + 1;
  int bin_q3 = i % nbins_q3 + 1;
  if (bin_q0 > nbins_q0) continue;

  double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
  double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
  if (q0 > q3) continue;

  h_param_frozen->SetBinContent(bin_q3, bin_q0, xsec->getParameterValue(i));
}


  std::vector<double> sum(nParams, 0.0);
  std::vector<double> sq_sum(nParams, 0.0);
  std::vector<double> val(nParams, 0.0);

  for (int i = 0; i < nParams; ++i) {
    post->SetBranchAddress(("xsec_" + std::to_string(i)).c_str(), &val[i]);
  }

  int nEntries = post->GetEntries();
  for (int e = 0; e < nEntries; ++e) {
    post->GetEntry(e);
    for (int i = 0; i < nParams; ++i) {
      sum[i] += val[i];
      sq_sum[i] += val[i] * val[i];
    }
  }

  // === Fill histograms ===
  for (int i = 0; i < nParams; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;

    if (bin_q0 > nbins_q0 || bin_q3 > nbins_q3) continue;
    double q0 = h_param_mean->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_param_mean->GetXaxis()->GetBinCenter(bin_q3);
    if (q0 > q3) continue;  // skip unphysical

    double mean = sum[i] / nEntries;
    double mean_sq = sq_sum[i] / nEntries;
    double variance = mean_sq - mean * mean;
    double sigma = variance > 0 ? std::sqrt(variance) : 0.0;

    h_param_mean->SetBinContent(bin_q3, bin_q0, mean);
    h_param_stddev->SetBinContent(bin_q3, bin_q0, sigma);
  }

  // === Output ===
  TCanvas* c = new TCanvas("c", "", 900, 700);
  gStyle->SetOptStat(0);
  c->Print("param_postfit_summary.pdf[");

  h_eventrate->Draw("COLZ");
  c->Print("param_postfit_summary.pdf");

  h_param_initn->Draw("COLZ");
  c->Print("param_postfit_summary.pdf");

  h_param_frozen->Draw("COLZ");
  c->Print("param_postfit_summary.pdf");

  h_param_mean->Draw("COLZ");
  c->Print("param_postfit_summary.pdf");

  h_param_stddev->Draw("COLZ");
  c->Print("param_postfit_summary.pdf");

  c->Print("param_postfit_summary.pdf]");

  TFile* fout = new TFile("param_postfit_summary.root", "RECREATE");
  h_eventrate->Write();
  h_param_init->Write();
  h_param_frozen->Write();
  h_param_mean->Write();
  h_param_stddev->Write();
  fout->Close();

  std::cout << "Saved: param_postfit_summary.pdf and .root" << std::endl;
  return 0;
}
