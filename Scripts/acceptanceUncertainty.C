#include <iostream>
#include <vector>
#include <cmath>
#include "DUNEStyle.h"

void acceptanceUncertainty() {
  // Get root file
  std::string filename = "Outputs/projections_outputs/Projections_CC_BField0_3_FidRad_160_FidLen_209_InstRad_249_45_InstLen_259_WithECal_seccurve_2layers_2Mev_2perc.root";
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

  double n_years = 6.;

  std::vector<double> x, y, exl, exh, eyl, eyh;
  
  TH1D* frac_uncertainty = new TH1D("N_uncertainty", "Fractional uncertainty due to acceptance correction", TrueEnuHist->GetNbinsX(), TrueEnuHist->GetXaxis()->GetXbins()->GetArray());
  for (int enubin=1; enubin<=TrueEnuHist->GetNbinsX(); enubin++) {
    double n_raw = TrueEnuHist->GetBinContent(enubin);
    double n_acc = AcceptedEnuHist->GetBinContent(enubin);
    double n_data = n_years*n_acc; // Assuming current simulations account for one year of data
    double acceptance = n_acc/n_raw;
    double acceptance_error = 0.01; // Dummy error of 1% on all acceptance bins due to model differences

    double N_min = std::round(n_data/10);
    double N_max = std::round(n_data*10);

    acceptance += acceptance_error;
    double N_lower = (1+2*n_data-acceptance-std::sqrt((1+2*n_data-acceptance)*(1+2*n_data-acceptance)-4*n_data*n_data))/(2*acceptance);
    acceptance -= 2*acceptance_error;
    double N_upper = (1+2*n_data-acceptance+std::sqrt((1+2*n_data-acceptance)*(1+2*n_data-acceptance)-4*n_data*n_data))/(2*acceptance);
    acceptance += acceptance_error;

    std::cout << "N_lower: " << N_lower << ", N_hat: " << n_data/acceptance << ", N_upper: " << N_upper << std::endl;

    frac_uncertainty->SetBinContent(enubin, (N_upper - N_lower)/n_data);
    
    x.push_back(TrueEnuHist->GetBinCenter(enubin));
    y.push_back(n_raw*n_years);
    double bw = TrueEnuHist->GetBinWidth(enubin)/2.;
    // exl.push_back(bw);
    // exh.push_back(bw);
    exl.push_back(0);
    exh.push_back(0);
    eyl.push_back(n_raw*n_years - N_lower);
    eyh.push_back(N_upper - n_raw*n_years);
  }
  
  // Create TCanvas
  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas("canvas", "Acceptance Correction Plots", 800, 600);
  gPad->SetTopMargin(0.13);
  const char* outputfilename = "Acceptance_Uncertainty_Plots.pdf";
  canvas->Print(Form("%s[", outputfilename));

  // Draw raw Enu hist with errors
  TrueEnuHist->Scale(n_years);
  TrueEnuHist->Draw("HIST");
  TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(), x.data(), y.data(), exl.data(), exh.data(), eyl.data(), eyh.data());
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.2);
  g->Draw("P same");
  canvas->Print(outputfilename);

  // Draw hist of fractional errors
  frac_uncertainty->Draw("HIST");
  canvas->Print(outputfilename);

  canvas->Print(Form("%s]", outputfilename));
  delete canvas;
  delete frac_uncertainty;
  delete g;
  rootfile->Close();
  delete rootfile;
}
