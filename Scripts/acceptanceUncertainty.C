#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Math/DistFunc.h"
#include "DUNEStyle.h"

enum class TestSide { TwoSided, Upper, Lower };

// Check whether observed n is consistent with N at confidence level (1 - alpha)
bool isConsistent(int n, int N, double a, double alpha, TestSide side) {
  double mu = N * a;
  double sigma = std::sqrt(N * a * (1 - a));
  // continuity correction
  double z = (n + 0.5 - mu) / sigma;
  double cdf = 0.5 * (1 + TMath::Erf(z / std::sqrt(2)));

  switch (side) {
    case TestSide::TwoSided:
      {
        double cdf_low  = cdf;            // P(n <= n_obs)
        double cdf_high = 1.0 - cdf_low;  // P(n >= n_obs)
        return (cdf_low >= alpha/2) && (cdf_high >= alpha/2);
      }
    case TestSide::Upper:
      return cdf >= alpha;  // only check lower tail: n_obs not too small
    case TestSide::Lower:
      return (1.0 - cdf) >= alpha; // only check upper tail: n_obs not too large
  }
  return false;
}

void acceptanceUncertainty() {
  // Get root file
  std::string filename = "Outputs/projections_outputs/Projections_CC_BField0_5_FidRad_160_FidLen_209_InstRad_249_45_InstLen_259_WithECal_seccurve_2layers_2Mev_2perc.root";
  TFile* rootfile = TFile::Open(filename.c_str(), "READ");
  if (!rootfile || rootfile->IsZombie()) {
    std::cerr << "Failed to open " << filename << std::endl;
  }

  // Get hists from root file
  TH1D* TrueEnuHist = dynamic_cast<TH1D*>(rootfile->Get("FHC_numu_NDGAr_TrueNeutrinoEnergy"));
  TH1D* AcceptedEnuHist = dynamic_cast<TH1D*>(rootfile->Get("FHC_numu_NDGAr_Accepted_TrueNeutrinoEnergy"));
  if (!TrueEnuHist || !AcceptedEnuHist) {
    std::cerr << "Could not find TrueEnuHist in " << filename << std::endl;
  }

  TH1D* AcceptanceHist = (TH1D*)AcceptedEnuHist->Clone("AcceptanceHist");
  AcceptanceHist->Divide(TrueEnuHist);
  AcceptanceHist->SetMaximum(1);

  int n_models = 3;
  int n_enu_bins = TrueEnuHist->GetNbinsX();
  std::vector<std::vector<double>> acceptance(n_models, std::vector<double>(n_enu_bins)); // Acceptance error in each bin

  double firstbinerr = 0.005;
  double lastbinerr = 0.02;
  for (int enubin = 1; enubin <= n_enu_bins; enubin++) {
    double acceptance_error = firstbinerr + (lastbinerr-firstbinerr)*(enubin-1)/(n_enu_bins-1);
    acceptance[0][enubin-1] = std::min(1., AcceptanceHist->GetBinContent(enubin) - acceptance_error);
    acceptance[1][enubin-1] = AcceptanceHist->GetBinContent(enubin);
    acceptance[2][enubin-1] = std::max(0., AcceptanceHist->GetBinContent(enubin) + acceptance_error);
    std::cout << enubin << ") acceptance = " << acceptance[1][enubin-1] << std::endl;
  }

  double n_years = 6.;

  std::vector<TH1D*> N_hists(n_models);
  std::vector<TH1D*> a_hists(n_models);
  std::vector<TH1D*> uncertainty_hists(n_models);
  std::vector<TGraphAsymmErrors*> error_bars(n_models);

  std::vector<double> maxN(n_enu_bins, 0.);
  std::vector<double> minN(n_enu_bins, 999999999.);
  std::vector<double> avgN(n_enu_bins, 0.);

  for (int model=0; model<n_models; model++) {
    std::vector<double> x, y, exl, exh, eyl, eyh;
    N_hists[model] = (TH1D*)TrueEnuHist->Clone("trueenuhist_model");
    a_hists[model] = (TH1D*)AcceptanceHist->Clone("acceptancehist_model");
    uncertainty_hists[model] = (TH1D*)TrueEnuHist->Clone("uncertaintyhist_model");

    for (int enubin=1; enubin<=n_enu_bins; enubin++) {
      double n_acc = AcceptedEnuHist->GetBinContent(enubin);
      double n_data = n_years*n_acc; // Assuming current simulations account for one year of data
      double a = acceptance[model][enubin-1];
      double N_hat = n_data/a;

      double N_lower = (1+2*n_data-a - std::sqrt((1+2*n_data-a)*(1+2*n_data-a) - 4*n_data*n_data))/(2*a);
      double N_upper = (1+2*n_data-a + std::sqrt((1+2*n_data-a)*(1+2*n_data-a) - 4*n_data*n_data))/(2*a);
      std::cout << "N_lower: " << N_lower << ", N_hat: " << N_hat << ", N_upper: " << N_upper << std::endl;

      N_hists[model]->SetBinContent(enubin, N_hat);
      a_hists[model]->SetBinContent(enubin, a);
      uncertainty_hists[model]->SetBinContent(enubin, (N_upper - N_lower)/N_hat);

      x.push_back(TrueEnuHist->GetBinCenter(enubin));
      y.push_back(N_hat);
      double bw = TrueEnuHist->GetBinWidth(enubin)/2.;
      exl.push_back(0);
      exh.push_back(0);
      eyl.push_back(N_hat - N_lower);
      eyh.push_back(N_upper - N_hat);

      if (N_upper > maxN[enubin-1]) maxN[enubin-1] = N_upper;
      if (N_lower < minN[enubin-1]) minN[enubin-1] = N_lower;
      avgN[enubin-1] += N_hat;
    }
    error_bars[model] = new TGraphAsymmErrors(x.size(), x.data(), y.data(), exl.data(), exh.data(), eyl.data(), eyh.data());
  }

  // Create TCanvas
  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas("canvas", "Acceptance Correction Plots", 800, 600);
  gPad->SetTopMargin(0.13);
  const char* outputfilename = "Acceptance_Uncertainty_Plots.pdf";
  canvas->Print(Form("%s[", outputfilename));
  int colourscheme[3] = {kOrange, kCyan-2, kRed+2};

  // Draw accepted Enu hist
  AcceptedEnuHist->Scale(n_years);
  AcceptedEnuHist->SetTitle("NDGAr Accepted Events");
  AcceptedEnuHist->Draw("HIST");
  canvas->Print(outputfilename);

  // Draw acceptance hist
  auto acceptance_leg = new TLegend(0.2, 0.2, 0.34, 0.35);
  for (int model=0; model<n_models; model++) {
    a_hists[model]->SetLineColor(colourscheme[model]);
    acceptance_leg->AddEntry(a_hists[model], ("Model " + std::to_string(model+1)).c_str(), "l");
    if (model==0) {
      a_hists[model]->SetTitle("NDGAr Acceptance");
      a_hists[model]->Draw("HIST");
    }
    else a_hists[model]->Draw("HIST SAME");
  }
  acceptance_leg->Draw();
  canvas->Print(outputfilename);

  // Draw total Enu hist with errors
  auto tot_events_leg = new TLegend(0.68, 0.67, 0.82, 0.82);
  for (int model=0; model<n_models; model++) {
    N_hists[model]->SetLineColor(colourscheme[model]);
    tot_events_leg->AddEntry(N_hists[model], ("Model " + std::to_string(model+1)).c_str(), "l");
    if (model==0) {
      N_hists[model]->SetTitle("NDGAr Total Events");
      N_hists[model]->Draw("HIST");
    }
    else N_hists[model]->Draw("HIST SAME");

    error_bars[model]->SetLineColor(kBlack);
    error_bars[model]->SetMarkerSize(0.2);
    // error_bars[model]->Draw("P same");
  }
  tot_events_leg->Draw();
  canvas->Print(outputfilename);

  // Draw hist of statistical fractional errors
  auto uncertainty_leg = new TLegend(0.2, 0.67, 0.34, 0.82);
  for (int model=0; model<n_models; model++) {
    uncertainty_hists[model]->SetLineColor(colourscheme[model]);
    uncertainty_leg->AddEntry(uncertainty_hists[model], ("Model " + std::to_string(model+1)).c_str(), "l");
    if (model==0) {
      uncertainty_hists[model]->SetTitle("Fractional binomial uncertainty on N due to acceptance correction");
      uncertainty_hists[model]->GetYaxis()->SetTitle("Fractional uncertainty");
      uncertainty_hists[model]->Draw("HIST");
    }
    else uncertainty_hists[model]->Draw("HIST SAME");
  }
  uncertainty_leg->Draw();
  canvas->Print(outputfilename);

  // Draw hist of model-difference fractional errors
  TH1D* modelErrorHist = new TH1D("modelErrorHist", "Fractional uncertainty on N from model differences; TrueNeutrinoEnergy; Fractional Uncertainty", n_enu_bins, TrueEnuHist->GetXaxis()->GetXbins()->GetArray());
  for (int enubin=1; enubin<=n_enu_bins; enubin++) {
    avgN[enubin-1] = avgN[enubin-1]/n_models;
    modelErrorHist->SetBinContent(enubin, (maxN[enubin-1] - minN[enubin-1]) / avgN[enubin-1]);
  }
  modelErrorHist->GetXaxis()->SetTitleSize(0.035);
  modelErrorHist->GetYaxis()->SetTitleSize(0.035);
  modelErrorHist->GetXaxis()->SetTitleOffset(1.2);
  modelErrorHist->GetYaxis()->SetTitleOffset(1.6);
  modelErrorHist->SetMaximum(0.1);
  modelErrorHist->SetLineColor(kGray+2);
  modelErrorHist->SetLineWidth(1);
  modelErrorHist->Draw("HIST");

  const int nPoints = 3;
  double x[nPoints] = {0.0, 2.5, 5.0};
  double y[nPoints] = {0.03, 0.03, 0.08}; // flat then linear rise
  
  auto stat_legend = new TLegend(0.2, 0.7, 0.54, 0.82);
  TGraph* line = new TGraph(nPoints, x, y);
  line->SetLineColor(kRed+1);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  stat_legend->AddEntry(line, "Phase II stat uncertainty", "l");
  line->Draw("L SAME");
  stat_legend->Draw();
  canvas->Print(outputfilename);

  canvas->Print(Form("%s]", outputfilename));
  delete canvas;
  rootfile->Close();
  delete rootfile;
}
