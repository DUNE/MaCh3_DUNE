void calc_roughness() {
  // Define the Enu bin edges and labels
  std::vector<std::pair<double, double>> enuBins = {
    {0.0, 1.5}, {1.5, 1.75}, {1.75, 2.0}, {2.0, 2.25}, {2.25, 2.5},
    {2.5, 2.75}, {2.75, 3.0}, {3.0, 3.25}, {3.25, 3.5}, {3.5, 4.0},
    {4.0, 4.5}, {4.5, 5.0}
  };

  // Create histogram to store roughness vs Enu
  TH1D* h_roughness_vs_enu = new TH1D("h_roughness_vs_enu", 
      "Roughness vs Enu;TrueNeutrinoEnergy (GeV);Roughness", 
      enuBins.size(), 0, enuBins.size());

  for (size_t i = 0; i < enuBins.size(); ++i) {
    auto [lo, hi] = enuBins[i];
    std::string filename = Form("Outputs/projections_outputs/loop/Projections_CC_Enu_%.2fto%.2f.root", lo, hi);
    TFile* fin = TFile::Open(filename.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
      std::cerr << "Failed to open " << filename << "\n";
      continue;
    }

    TH2D* h_raw = dynamic_cast<TH2D*>(fin->Get("FHC_numu_NDGAr_TrueQ3_vs_TrueQ0"));
    TH2D* h_accepted = dynamic_cast<TH2D*>(fin->Get("FHC_numu_NDGAr_Accepted_TrueQ3_vs_TrueQ0"));
    if (!h_raw || !h_accepted) {
      std::cerr << "Histograms missing in " << filename << "\n";
      fin->Close();
      continue;
    }

    // Clone to divide safely
    TH2D* h_acceptance = dynamic_cast<TH2D*>(h_accepted->Clone("h_acceptance"));
    h_acceptance->Divide(h_raw);
    
    for (int x = 1; x <= h_raw->GetNbinsX(); ++x) {
      for (int y = 1; y <= h_raw->GetNbinsY(); ++y) {
        int binnum = h_acceptance->GetBin(x, y);
        double num = h_accepted->GetBinContent(x, y);
        double den = h_raw->GetBinContent(x, y);
        if (den == 0) h_acceptance->SetBinContent(binnum, 0);
        else if (num == 0) h_acceptance->SetBinContent(binnum, 1e-6);
      }
    }

    // Compute roughness
    double roughness = 0;
    for (int x = 1; x <= h_acceptance->GetNbinsX(); ++x) {
      for (int y = 1; y <= h_acceptance->GetNbinsY(); ++y) {
        double centerVal = h_acceptance->GetBinContent(x, y);
        if (centerVal <= 0) continue;

        for (int dx = -1; dx <= 1; ++dx) {
          for (int dy = -1; dy <= 1; ++dy) {
            if (dx == 0 && dy == 0) continue;
            int nx = x + dx;
            int ny = y + dy;
            if (nx < 1 || nx > h_acceptance->GetNbinsX()) continue;
            if (ny < 1 || ny > h_acceptance->GetNbinsY()) continue;

            double neighborVal = h_acceptance->GetBinContent(nx, ny);
            if (neighborVal <= 0) continue;

            double diff = centerVal - neighborVal;
            roughness += diff * diff;
          }
        }
      }
    }

    h_roughness_vs_enu->SetBinContent(i + 1, roughness);
    h_roughness_vs_enu->GetXaxis()->SetBinLabel(i + 1, Form("%.2f–%.2f", lo, hi));
    fin->Close();
  }

  // Write output to file
  TFile* fout = new TFile("roughness_vs_enu.root", "RECREATE");
  h_roughness_vs_enu->Write();
  fout->Close();
}

