#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <TLatex.h>

void drawPage(TCanvas* c, TH2D* h_nom, TH2D* h_var, TH2D* h_ratio,
              const std::string& title, const std::string& varLabel)
{
    c->Clear();
    c->Divide(3,1);

    c->cd(0);

    TLatex latex;
    latex.SetTextFont(62);   // bold Helvetica-like
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.5, 0.96, title.c_str());

    // centered at top of canvas

    

    // ---- Nominal ----
    c->cd(1);
    gPad->SetRightMargin(0.15);
    h_nom->SetTitle((title + " Nominal").c_str());
    h_nom->Draw("colz");

    // ---- Variation ----
    c->cd(2);
    gPad->SetRightMargin(0.15);
    h_var->SetTitle((title + " " + varLabel).c_str());
    h_var->Draw("colz");

    // ---- Ratio ----
    c->cd(3);
    gPad->SetRightMargin(0.15);
   h_ratio->SetContour(100);
    h_ratio->SetTitle((varLabel + " fractional variation").c_str());
    h_ratio->GetZaxis()->SetTitle("(Var - Nom) / Nom");
    h_ratio->GetYaxis()->SetTitle("Reconstructed Neutrino Energy [GeV]");
    //h_ratio->GetYaxis()->SetRangeUser(0,5.0);
    //h_ratio->GetXaxis()->SetRangeUser(0,5.0);
    h_ratio->GetXaxis()->SetTitle("True Neutrino Energy [GeV]");
    h_ratio->SetMinimum(-0.3);
    h_ratio->SetMaximum(+0.3);
    h_ratio->Draw("colz");

    
}

void makeSigmaPDF_2D() {

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Smooth gradients
    gStyle->SetNumberContours(100);

    TFile* f = TFile::Open("/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/sigmavar_recovtrue.root");
    //TFile* f = TFile::Open("/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/SigmaVariation_150426_template_recobinning.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error opening file\n";
        return;
    }

    const Int_t NRGBs = 3;
    const Int_t NCont = 100;

    Double_t stops[NRGBs] = {0.0, 0.5, 1.0};

    // Blue → White → Red
    Double_t red[NRGBs]   = {0.0, 1.0, 1.0};
    Double_t green[NRGBs] = {0.0, 1.0, 0.0};
    Double_t blue[NRGBs]  = {1.0, 1.0, 0.0};

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    TCanvas* c = new TCanvas("c","c",1200,400);

    std::string pdfName = "SigmaVariations_2D_split_template_truth_friday.pdf";
    c->Print((pdfName + "[").c_str()); // open PDF

    TIter next(f->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {

        if (std::string(key->GetClassName()) != "TDirectoryFile") continue;

        TDirectory* systDir = (TDirectory*)key->ReadObj();
        std::string systName = systDir->GetName();

        std::cout << "Processing: " << systName << std::endl;

        // Take first channel (same assumption as before)
        TIter next2(systDir->GetListOfKeys());
        TKey* key2 = (TKey*)next2();
        if (!key2) continue;

        TDirectory* chanDir = (TDirectory*)key2->ReadObj();
        std::string chanName = chanDir->GetName();

        // Get histograms
        TH2D* h_nom = (TH2D*)chanDir->Get("Variation_2");
        TH2D* h_dn  = (TH2D*)chanDir->Get("Variation_1");
        TH2D* h_up  = (TH2D*)chanDir->Get("Variation_3");

        if (!h_nom || !h_dn || !h_up) {
            std::cout << "  Missing histograms\n";
            continue;
        }

        // Clone ratios
        TH2D* r_up = (TH2D*)h_up->Clone((systName+"_fup").c_str());
        TH2D* r_dn = (TH2D*)h_dn->Clone((systName+"_fdn").c_str());

        // (var - nom)/nom
        r_up->Add(h_nom, -1.0);
        r_dn->Add(h_nom, -1.0);

        for (int i = 1; i <= h_nom->GetNbinsX(); ++i) {
            for (int j = 1; j <= h_nom->GetNbinsY(); ++j) {
                if (h_nom->GetBinContent(i,j) == 0) {
                    r_up->SetBinContent(i,j, 0);
                    r_dn->SetBinContent(i,j, 0);
                }
            }
        }

        r_up->Divide(h_nom);
        r_dn->Divide(h_nom);



        std::string title = systName + " (" + chanName + ")";

        // ---- PAGE 1: +1σ ----
        drawPage(c, h_nom, h_up, r_up, title, "+1#sigma");
        c->Print(pdfName.c_str());

        // ---- PAGE 2: -1σ ----
        drawPage(c, h_nom, h_dn, r_dn, title, "-1#sigma");
        c->Print(pdfName.c_str());
    }

    c->Print((pdfName + "]").c_str()); // close PDF

    std::cout << "Saved to " << pdfName << std::endl;
}