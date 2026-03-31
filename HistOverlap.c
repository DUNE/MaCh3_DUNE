#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"


void overlay() {
    // Open files
    //TFile* OscMF = TFile::Open("OscTrueParamsOnRecoNoNC.root");
    TFile* OscMF = TFile::Open("OscMFNoNC.root");
    TFile* OscPMNS = TFile::Open("OscPMNSNoNC.root");

    TIter next(OscMF->GetListOfKeys()); 
    TKey* key; 
    while ((key = (TKey*)next())) { 
        auto HistoOsc = OscMF->Get<TH1D>(key->GetName());
        auto HistoUnosc = OscPMNS->Get<TH1D>(key->GetName());

        // Style
        HistoOsc->SetLineColor(kRed);
        HistoOsc->SetLineWidth(2);

        HistoUnosc->SetLineColor(kBlue);
        HistoUnosc->SetLineWidth(2);

        HistoOsc->SetStats(0);
        HistoUnosc->SetStats(0);

        // Draw
        TCanvas *c = new TCanvas("c","Overlay");
        HistoUnosc->Draw("HIST");
        HistoOsc->Draw("HIST SAME");

        // Add legend
        auto leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->AddEntry(HistoOsc,"Model Free","l");
        leg->AddEntry(HistoUnosc,"PMNS","l");
        leg->Draw();

        c->SaveAs(Form("Histo_%s.png", key->GetName()));
    }
}