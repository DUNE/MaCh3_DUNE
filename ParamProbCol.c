#include <vector>
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

const int NParams = 1920;

std::vector<double> BinStart = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0,6.25,6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25,8.5,8.75,9.0,9.25,9.5,9.75,10.0};
std::vector<double> BinMid = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,6.125,6.375,6.625,6.875,7.125,7.375,7.625,7.875,8.125,8.375,8.625,8.875,9.125,9.375,9.625,9.875};

std::vector<double> NormParams(NParams, 0.0);

void NormaliseAndPlot(TH2D* h, TString HistName)
{
    for (int x = 1; x <= h->GetNbinsX(); x++) {
        double ColSum = 0;
        for (int y = 1; y <= h->GetNbinsY(); y++) {
            ColSum += h->GetBinContent(x,y);
        }
        if (ColSum == 0) continue;
        for (int y = 1; y <= h->GetNbinsY(); y++) {
            double Content = h->GetBinContent(x,y);
            h->SetBinContent(x, y, Content/ColSum);
        }
    }

    auto c = new TCanvas();
    c->Draw();
    c->cd();
    h->Draw();
    h->GetXaxis()->SetTitle("Energy (GeV)");
    h->GetYaxis()->SetTitle("Oscillation Probability");
    h->Draw("COLZ");
    h->GetZaxis()->SetTitle("Number of Events");
    h->GetZaxis()->SetTitleOffset(1.5);
    h->SetStats(0);
    c->SetRightMargin(0.15);
    c->Update();
    c->Print(HistName + "OscProbCol.pdf");
}

void treeApp_read(TString file_name)
{
    auto f = TFile::Open(file_name);

    if(f->IsZombie() || !f){
        std::cerr<<"Error cannot find "<<file_name<<std::endl;
        throw;
    }

    auto t1 = f->Get<TTree>("posteriors");

    for(int i = 0; i < NParams; i++){
        std::string name = "xsec_" + std::to_string(i); 
        t1->SetBranchAddress(name.c_str(), &NormParams[i]);
    }

    auto HistFHCnumu = new TH2D("HistFHCnumu", "Oscillation Probability for FHC_numu Sample", NParams/48, BinStart.data(), 20, 0.0, 1.0); 
    auto HistFHCnue = new TH2D("HistFHCnue", "Oscillation Probability for FHC_nue Sample", NParams/48, BinStart.data(), 20, 0.0, 1.0); 
    auto HistRHCnumu = new TH2D("HistRHCnumu", "Oscillation Probability for RHC_numu Sample", NParams/48, BinStart.data(), 20, 0.0, 1.0); 
    auto HistRHCnue = new TH2D("HistRHCnue", "Oscillation Probability for RHC_nue Sample", NParams/48, BinStart.data(), 20, 0.0, 1.0); 

    Long64_t nentries = t1->GetEntries();
    for (Long64_t i=2000000; i<nentries; i++) {
        t1->GetEntry(i);
        for(int p = 0; p < NParams; p++){
            int k = p % (NParams/48);
            if( p < NParams/4 ){
                HistFHCnumu->Fill(BinMid[k], NormParams[p]);
            }
            else if( p >= NParams/4 && p < NParams/2){
                HistFHCnue->Fill(BinMid[k], NormParams[p]);
            }
            else if( p >= NParams/2 && p < 3*NParams/4){
                HistRHCnumu->Fill(BinMid[k], NormParams[p]);
            }
            else{
                HistRHCnue->Fill(BinMid[k], NormParams[p]);
            }
        }
    }

    NormaliseAndPlot(HistFHCnumu, "FHCnumu");
    NormaliseAndPlot(HistFHCnue, "FHCnue");
    NormaliseAndPlot(HistRHCnumu, "RHCnumu");
    NormaliseAndPlot(HistRHCnue, "RHCnue");

    // for (int x = 1; x <= HistFHCnumu->GetNbinsX(); x++) {
    //     double ColSum = 0;
    //     for (int y = 1; y <= HistFHCnumu->GetNbinsY(); y++) {
    //         double Content = HistFHCnumu->GetBinContent(x,y);
    //         ColSum += Content;
    //     }
    //     if (ColSum == 0) continue;
    //     for (int y = 1; y <= HistFHCnumu->GetNbinsY(); y++) {
    //         double Content = HistFHCnumu->GetBinContent(x,y);
    //         HistFHCnumu->SetBinContent(x, y, Content/ColSum);
    //     }
    // }

    // auto c = new TCanvas();
    // c->Draw();
    // c->cd();
    // HistFHCnumu->Draw();
    // HistFHCnumu->GetXaxis()->SetTitle("Energy (GeV)");
    // HistFHCnumu->GetYaxis()->SetTitle("Oscillation Probability");
    // HistFHCnumu->Draw("COLZ");
    // HistFHCnumu->GetZaxis()->SetTitle("Number of Events");
    // HistFHCnumu->GetZaxis()->SetTitleOffset(1.5);
    // HistFHCnumu->SetStats(0);
    // c->SetRightMargin(0.15);
    // c->Update();
    // c->Print("FHCnumuOscProbCol.pdf");
}