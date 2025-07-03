void movePalette(TH1* hist) {
  gPad->Update(); // palette doesn't exist until after first draw
  TPaletteAxis* palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
  if (palette) {
    palette->SetX1NDC(0.855);  // Left edge
    palette->SetX2NDC(0.9);  // Right edge
  }
}

void changeAxisTitle(TAxis* axis) {
  std::string title = axis->GetTitle();
  if (title ==  "TrueQ0") axis->SetTitle("q_{0} [GeV]");
  else if (title == "TrueQ3") {
    axis->SetTitle("q_{3} [GeV/c]");
    axis->SetNdivisions(4,8,0);
  }
}

double getMaximumMagnitudeBinVal(TH1* hist) {
  double maxval = hist->GetMaximum();
  double minval = hist->GetMinimum();

  if (std::abs(minval) > maxval) return std::abs(minval);
  else return maxval;
}

void geom_correction_plots() {
  // Open the ROOT files
  TFile* f_nom = TFile::Open("outputs/Projections_CC_BField0_5_FidRad_160_FidLen_209_InstRad_249_45_InstLen_259_WithECal_Enu_1to5GeV_clean.root");
  TFile* f_geom = TFile::Open("outputs/Projections_CC_BField0_5_FidRad_160_FidLen_209_InstRad_249_45_InstLen_259_WithECal_Enu_1to5GeV_geomcorr_clean.root");

  if (!f_nom || f_nom->IsZombie() || !f_geom || f_geom->IsZombie()) {
    std::cerr << "Error: Failed to open one or both ROOT files." << std::endl;
    return;
  }

  // Load histograms
  TH2D* h1 = dynamic_cast<TH2D*>(f_nom->Get("FHC_numu_NDGAr_TrueQ3_vs_TrueQ0"));
  TH2D* h2 = dynamic_cast<TH2D*>(f_nom->Get("FHC_numu_NDGAr_Accepted_TrueQ3_vs_TrueQ0"));
  TH2D* g1 = dynamic_cast<TH2D*>(f_geom->Get("FHC_numu_NDGAr_TrueQ3_vs_TrueQ0"));
  TH2D* g2 = dynamic_cast<TH2D*>(f_geom->Get("FHC_numu_NDGAr_Accepted_TrueQ3_vs_TrueQ0"));

  if (!h1 || !h2 || !g1 || !g2) {
    std::cerr << "Error: Could not find all required histograms." << std::endl;
    return;
  }

  // Rebin each histogram by 2x2
  h1->Rebin2D(2, 2);
  h2->Rebin2D(2, 2);
  g1->Rebin2D(2, 2);
  g2->Rebin2D(2, 2);

  // Compute efficiency histograms
  TH2D* a1 = dynamic_cast<TH2D*>(h2->Clone("a1"));
  a1->Divide(h1);

  TH2D* a2 = dynamic_cast<TH2D*>(g2->Clone("a2"));
  a2->Divide(g1);

  // Compute difference a2 - a1
  TH2D* a2_over_a1 = dynamic_cast<TH2D*>(a2->Clone("a2_over_a1"));
  a2_over_a1->Add(a1, -1);
  int nx = a2_over_a1->GetNbinsX();
  int ny = a2_over_a1->GetNbinsY();

  // Compute (g1 - h1) / sqrt(g1 + h1)
  TH2D* diff = dynamic_cast<TH2D*>(g1->Clone("diff"));
  diff->Add(h1, -1); // g1 - h1

  nx = diff->GetNbinsX();
  ny = diff->GetNbinsY();

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      double gval = g1->GetBinContent(ix, iy);
      double hval = h1->GetBinContent(ix, iy);
      double denom = std::sqrt(2*hval);
      if (denom > 0) {
        diff->SetBinContent(ix, iy, (gval - hval) / denom);
      } else {
        diff->SetBinContent(ix, iy, 0);
      }
    }
  }

  // Plot everything to a multipage PDF
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("1.2f");
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);

  const Int_t NRGBs = 3;
  const Int_t NCont = 255;

  // Define positions (0 to 1) and RGB colors
  Double_t stops[NRGBs]    = { 0.00, 0.5, 1.00 };
  Double_t red[NRGBs]      = { 0.00, 1.0, 1.00 };  // Blue -> White -> Red
  Double_t green[NRGBs]    = { 0.00, 1.0, 0.00 };
  Double_t blue[NRGBs]     = { 1.00, 1.0, 0.00 };

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);


  TCanvas* c = new TCanvas("c", "", 900, 600);
  c->Print("efficiency_comparison_TH2D.pdf[");

  c->cd(); a1->SetTitle("Acceptance (Nominal: h2 / h1)");
  changeAxisTitle(a1->GetXaxis());
  changeAxisTitle(a1->GetYaxis());
  a1->Draw("COLZ");
  a1->GetXaxis()->SetLabelSize(0.04);
  a1->GetXaxis()->SetLabelFont(42);
  a1->GetXaxis()->SetTitleOffset(1.3);
  a1->GetYaxis()->SetLabelSize(0.04);
  a1->GetYaxis()->SetLabelFont(42);
  gPad->SetFrameLineWidth(2);
  movePalette(a1);
  c->Print("efficiency_comparison_TH2D.pdf");

  c->cd(); a2->SetTitle("Acceptance (Geom-Corrected: g2 / g1)");
  changeAxisTitle(a2->GetXaxis());
  changeAxisTitle(a2->GetYaxis());
  a2->Draw("COLZ");
  a2->GetXaxis()->SetLabelSize(0.04);
  a2->GetXaxis()->SetLabelFont(42);
  a2->GetXaxis()->SetTitleOffset(1.3);
  a2->GetYaxis()->SetLabelSize(0.04);
  a2->GetYaxis()->SetLabelFont(42);
  gPad->SetFrameLineWidth(2);
  movePalette(a2);
  c->Print("efficiency_comparison_TH2D.pdf");

  c->cd(); a2_over_a1->SetTitle("Acceptance Difference: a2 - a1");
  changeAxisTitle(a2_over_a1->GetXaxis());
  changeAxisTitle(a2_over_a1->GetYaxis());
  a2_over_a1->Draw("COLZ");
  a2_over_a1->GetXaxis()->SetLabelSize(0.04);
  a2_over_a1->GetXaxis()->SetLabelFont(42);
  a2_over_a1->GetXaxis()->SetTitleOffset(1.3);
  a2_over_a1->GetYaxis()->SetLabelSize(0.04);
  a2_over_a1->GetYaxis()->SetLabelFont(42);
  double max = getMaximumMagnitudeBinVal(a2_over_a1);
  a2_over_a1->SetMinimum(-max);
  a2_over_a1->SetMaximum(max);
  gPad->SetFrameLineWidth(2);
  movePalette(a2_over_a1);
  c->Print("efficiency_comparison_TH2D.pdf");

  c->cd(); diff->SetTitle("Event rates difference: (g1 - h1) / sqrt(2g1)");
  changeAxisTitle(diff->GetXaxis());
  changeAxisTitle(diff->GetYaxis());
  diff->Draw("COLZ"); 
  diff->GetXaxis()->SetLabelSize(0.04);
  diff->GetXaxis()->SetLabelFont(42);
  diff->GetXaxis()->SetTitleOffset(1.3);
  diff->GetYaxis()->SetLabelSize(0.04);
  diff->GetYaxis()->SetLabelFont(42);
  max = getMaximumMagnitudeBinVal(diff);
  diff->SetMinimum(-max);
  diff->SetMaximum(max);
  gPad->SetFrameLineWidth(2);
  movePalette(diff);
  c->Print("efficiency_comparison_TH2D.pdf");

  c->Print("efficiency_comparison_TH2D.pdf]");
}

