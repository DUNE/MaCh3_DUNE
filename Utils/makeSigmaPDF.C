#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

TH1D* getProjection(TDirectory* dir, const char* name, const std::string& tag) {
    TH2D* h2 = (TH2D*)dir->Get(name);
    if (!h2) return nullptr;
    return h2->ProjectionX(tag.c_str());
}

void makeSigmaPDF() {
    TFile* f = TFile::Open("/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/sigmavar_recovtrue.root");

    if (!f || f->IsZombie()) {
        std::cout << "Error opening file\n";
        return;
    }

    TCanvas* c = new TCanvas("c","c",800,800);

    std::string pdfName = "SigmaVariations_template.pdf";
    c->Print((pdfName + "[").c_str()); // open PDF

    TIter next(f->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {

        if (std::string(key->GetClassName()) != "TDirectoryFile") continue;

        TDirectory* systDir = (TDirectory*)key->ReadObj();
        std::string systName = systDir->GetName();

        std::cout << "Processing systematic: " << systName << std::endl;

        // 👉 Take FIRST channel (you can change this later)
        TIter next2(systDir->GetListOfKeys());
        TKey* key2 = (TKey*)next2();
        if (!key2) continue;

        TDirectory* chanDir = (TDirectory*)key2->ReadObj();
        std::string chanName = chanDir->GetName();

        // Get correct variations
        TH1D* h_nom = getProjection(chanDir, "Variation_2", systName+"_nom");
        TH1D* h_dn  = getProjection(chanDir, "Variation_1", systName+"_dn");
        TH1D* h_up  = getProjection(chanDir, "Variation_3", systName+"_up");

        if (!h_nom || !h_dn || !h_up) {
            std::cout << "  Missing histograms\n";
            continue;
        }

        // Ratios
        TH1D* r_up = (TH1D*)h_up->Clone((systName+"_rup").c_str());
        TH1D* r_dn = (TH1D*)h_dn->Clone((systName+"_rdn").c_str());

        r_up->Divide(h_nom);
        r_dn->Divide(h_nom);

        // Clear canvas
        c->Clear();

        // Pads
        TPad* p1 = new TPad("p1","",0,0.3,1,1);
        TPad* p2 = new TPad("p2","",0,0,1,0.3);

        p1->SetBottomMargin(0.02);
        p2->SetTopMargin(0.05);
        p2->SetBottomMargin(0.3);

        p1->Draw();
        p2->Draw();

        // ---- TOP PAD ----
        p1->cd();

        h_nom->SetLineColor(kBlack);
        h_up->SetLineColor(kRed);
        h_dn->SetLineColor(kBlue);

        h_nom->SetTitle((systName + " (" + chanName + ")").c_str());

        h_nom->Draw("hist");
        h_up->Draw("hist same");
        h_dn->Draw("hist same");

        TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->AddEntry(h_nom,"Nominal","l");
        leg->AddEntry(h_up,"+1#sigma","l");
        leg->AddEntry(h_dn,"-1#sigma","l");
        leg->Draw();

        // ---- RATIO ----
        p2->cd();

        r_up->SetLineColor(kRed);
        r_dn->SetLineColor(kBlue);

        r_up->SetMinimum(0.8);
        r_up->SetMaximum(1.2);

        r_up->SetTitle("");
        r_up->GetYaxis()->SetTitle("Ratio");

        r_up->Draw("hist");
        r_dn->Draw("hist same");

        TLine* line = new TLine(
            r_up->GetXaxis()->GetXmin(),1,
            r_up->GetXaxis()->GetXmax(),1
        );
        line->SetLineStyle(2);
        line->Draw();

        // Print page
        c->Print(pdfName.c_str());

        delete p1;
        delete p2;
    }

    c->Print((pdfName + "]").c_str()); // close PDF

    std::cout << "Saved to " << pdfName << std::endl;
}