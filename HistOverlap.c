#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <fmt/core.h>


void overlay() { // Create a script to print two histograms on same canvas
    TFile* OscMF = TFile::Open("PostPredPMNSv2.root"); // Load in histograms
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

        c->SaveAs("PostPred.png"); // Saving and creating title
    }
}


void PostPred() { // Create a script to print two histograms on same canvas
    TFile* OscPMNS = TFile::Open("PostPredPMNSv2.root"); // Load in histograms
    TFile* OscMF = TFile::Open("PostPredFHCnumuv2.root");
    TFile* OscData = TFile::Open("OscPMNSNoNC.root");

    auto HistoPMNS = OscPMNS->Get<TH1D>("Predictive/FD_FHC_numu/FD_FHC_numu_mc_PostPred"); // PostPredPMNSv2
    auto HistoMF = OscMF->Get<TH1D>("Predictive/FD_FHC_numu/FD_FHC_numu_mc_PostPred");
    auto HistoData = OscData->Get<TH1D>("hRecoNeutrinoEnergyFHC_numu");

    HistoPMNS->SetLineColor(kRed); // Cosmetics for PMNS model histogram
    HistoPMNS->SetLineWidth(2);
    HistoPMNS->SetStats(0);

    HistoMF->SetLineColor(kBlue); // Cosmetics for Model Free model histogram
    HistoMF->SetLineWidth(2);
    HistoMF->SetStats(0);

    HistoData->SetLineColor(kBlack); // Cosmetics for data histogram
    HistoData->SetLineWidth(2);
    HistoData->SetLineStyle(2);
    HistoData->SetStats(0);

    TCanvas *c = new TCanvas("c","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
    HistoMF->Draw("HIST E0");
    HistoMF->SetTitle("Posterior Predictive FHC_numu");
    HistoPMNS->Draw("HIST SAME E0");
    HistoData->Draw("HIST SAME");


    auto leg = new TLegend(0.5, 0.5, 0.88, 0.88); // Creating the legend
    leg->SetMargin(0.4);
    leg->AddEntry(HistoMF,fmt::format("Model Free - #int = {:.2f}", HistoMF->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoPMNS,fmt::format("PMNS - #int = {:.2f}", HistoPMNS->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoData,fmt::format("Data - #int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
    leg->Draw();

    c->SaveAs("PostPred.png"); // Saving and creating title

    // int NumBins = HistoModel->GetNbinsX(); // Find number of bins (energy normalisation parameters) for each histogram
    // for(int j = 1; j <= NumBins; j++) {
    //     double BinSizeModel = HistoModel->GetBinContent(j); // Find number of events for each energy parameter
    //     double BinSizeData = HistoData->GetBinContent(j);
    //     double Param;
    //     Param = BinSizeModel / BinSizeData;
    //     std::cout << Param << "\n";
    // }
}