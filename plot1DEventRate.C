void plot1DEventRate() {

  // Open file
  TFile *f = TFile::Open("/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/eventrate_reco.root");

  // Get 2D histogram
  TH2D *h2 = (TH2D*)f->Get("hOnAxis_numuCC_numode");

  // --------------------------------------------------
  // Make 1D event rate histogram
  // --------------------------------------------------
  // ProjectionX = integrate over Y
  // ProjectionY = integrate over X
  //
  // Since your X axis is:
  // Reconstructed Neutrino Energy [GeV]
  //
  // this gives event rate vs reconstructed neutrino energy
  // --------------------------------------------------

  TH1D *h1 = h2->ProjectionX("h1");

  // Optional: remove first bin
  //h1->GetXaxis()->SetRange(2, h1->GetNbinsX());
  h1->GetXaxis()->SetRangeUser(0, 10.0);

  // --------------------------------------------------
  // Style settings (publication style)
  // --------------------------------------------------

  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c","c",800,700);

  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);

  h1->SetTitle("");

  h1->GetXaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
  h1->GetYaxis()->SetTitle("Events");

  h1->GetXaxis()->SetTitleSize(0.045);
  h1->GetYaxis()->SetTitleSize(0.045);

  h1->GetXaxis()->SetLabelSize(0.040);
  h1->GetYaxis()->SetLabelSize(0.040);

  h1->GetYaxis()->SetTitleOffset(1.4);

  // Line style
  h1->SetLineWidth(3);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.8);

  // Draw with error bars
  h1->Draw("HIST E");

  // Optional log scale for event rates
  // gPad->SetLogy();

  // Save
  c->SaveAs("eventrate_1D.pdf");
}