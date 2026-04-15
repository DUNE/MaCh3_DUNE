#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

void overlay() { // Create a script to print two histograms on same canvas
    TFile* OscMF = TFile::Open("OscMFNoNC.root"); // Load in histograms
    TFile* OscPMNS = TFile::Open("OscPMNSNoNC.root");

    TIter next(OscMF->GetListOfKeys()); // Producing list of keys (histograms) from file
    TKey* key; 
    while ((key = (TKey*)next())) { 
        auto HistoOsc = OscMF->Get<TH1D>(key->GetName()); // Get titles of individual histograms
        auto HistoUnosc = OscPMNS->Get<TH1D>(key->GetName());

        HistoOsc->SetLineColor(kRed); // Cosmetics for oscillated histogram
        HistoOsc->SetLineWidth(2);
        HistoOsc->SetStats(0);

        HistoUnosc->SetLineColor(kBlue); // Cosmetics for unoscillated histogram
        HistoUnosc->SetLineWidth(2);
        HistoUnosc->SetStats(0);

        TCanvas *c = new TCanvas("c","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
        HistoUnosc->Draw("HIST");
        HistoOsc->Draw("HIST SAME");

        auto leg = new TLegend(0.65, 0.75, 0.88, 0.88); // Creating the legend
        leg->AddEntry(HistoOsc,"Model Free","l");
        leg->AddEntry(HistoUnosc,"PMNS","l");
        leg->Draw();

        c->SaveAs(Form("Histo_%s.png", key->GetName())); // Saving and creating title
    }
}