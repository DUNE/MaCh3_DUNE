void plotEventRate() {

  // Open file
  TFile *f = TFile::Open("/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/eventrate_reco.root");

  // Get histogram
  TH2D *h = (TH2D*)f->Get("hOnAxis_numuCC_numode");

  // Style settings
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(100);

  // Remove first bin by restricting axis range
  int firstBinX = 2;
  int lastBinX  = h->GetNbinsX();
  int firstBinY = 2;
  int lastBinY  = h->GetNbinsY();

  h->GetXaxis()->SetRange(firstBinX, lastBinX);
  h->GetYaxis()->SetRange(firstBinY, lastBinY);

  // Labels
  h->SetTitle("");
  // h->GetXaxis()->SetTitle("True Neutrino Energy [GeV]");
  // h->GetYaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
  h->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
  h->GetYaxis()->SetTitle("Reconstructed Lepton Energy [GeV]");
  h->GetZaxis()->SetTitle("Events");

  // Canvas
  TCanvas *c = new TCanvas("c","c",800,700);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);

  // Draw
  h->Draw("COLZ");

  // Save
  c->SaveAs("eventratereco.pdf");
}