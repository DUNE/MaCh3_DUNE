#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <fmt/core.h>
#include <unordered_map>
#include <string>
#include <iostream>

void overlay() { // Create a script to print two histograms on same canvas
    TFile* OscMF = TFile::Open("EventRates/OscPMNSNoNC.root"); // Load in histograms
    TFile* OscPMNS = TFile::Open("EventRates/UnoscNoNC.root");

    //TIter next(OscMF->GetListOfKeys()); // Producing list of keys (histograms) from file
    //TKey* key; 
    //while ((key = (TKey*)next())) { 
    auto HistoOsc = OscMF->Get<TH1D>("hRecoNeutrinoEnergyFHC_numu"); // Get titles of individual histograms
    auto HistoUnosc = OscPMNS->Get<TH1D>("hRecoNeutrinoEnergyFHC_numu"); // Get titles of individual histograms

    HistoOsc->SetLineColor(kBlack); // Cosmetics for oscillated histogram
    HistoOsc->SetLineWidth(2);
    HistoOsc->SetStats(0);

    HistoUnosc->SetLineColor(kMagenta); // Cosmetics for unoscillated histogram
    HistoUnosc->SetLineWidth(2);
    HistoUnosc->SetStats(0);

    TCanvas *c = new TCanvas("c","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
    HistoUnosc->Draw("HIST");
    HistoUnosc->SetTitle("Event Rate Comparison - FHC_numu");
    HistoUnosc->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    HistoUnosc->GetYaxis()->SetTitle("Number of Events");
    HistoOsc->Draw("HIST SAME");

    auto leg = new TLegend(0.6, 0.6, 0.88, 0.88); // Creating the legend
    leg->AddEntry(HistoUnosc,"FD Unoscillated","l");
    leg->AddEntry(HistoOsc,"FD Oscillated","l");
    leg->Draw();
    c->SaveAs("PostPred.png"); // Saving and creating title
    //}
}


void PostPred(TString PPFile) { // Create a script to print two histograms on same canvas
    std::string Sample = PPFile.Data();
    std::string PlotFile = Sample.substr(0, Sample.size() - 5);
    
    TFile* OscPMNS = TFile::Open("AllSampleFits/PostPredPMNS.root"); // Load in histograms
    TFile* OscMF = TFile::Open(PPFile);
    TFile* OscData = TFile::Open("EventRates/OscPMNSNoNC.root");

    auto HistoPMNS = OscPMNS->Get<TH1D>("Predictive/FHC_numu/FHC_numu_mc_PostPred"); // PostPredPMNSv2
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

    TCanvas *c = new TCanvas("PostPredDataAndMF","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
    HistoMF->Draw("HIST E0");
    HistoMF->SetTitle("Posterior Predictive - FHC_numu");
    //HistoPMNS->Draw("HIST SAME E0");
    //HistoPMNS->GetXaxis()->SetRangeUser(10, 30);
    HistoData->Draw("HIST SAME");
    HistoMF->GetXaxis()->SetRangeUser(0, 20);
    HistoMF->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    HistoMF->GetYaxis()->SetTitle("Number of Events");
    HistoData->GetXaxis()->SetRangeUser(0, 20);

    auto leg = new TLegend(0.5, 0.5, 0.88, 0.88); // Creating the legend
    leg->SetMargin(0.4);
    leg->AddEntry(HistoMF,fmt::format("Model Free - #int = {:.2f}", HistoMF->Integral()).c_str(),"lpf");
    //leg->AddEntry(HistoPMNS,fmt::format("PMNS - #int = {:.2f}", HistoPMNS->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoData,fmt::format("Data - #int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
    leg->Draw();

    //c->SaveAs((PlotFile + "Overlay.root").c_str());
    c->SaveAs((PlotFile + "(0,20).png").c_str()); // Saving and creating title
}

void RegValComp(TString Sam) { // Create a script to print two histograms on same canvas
    std::string Sample = Sam.Data();
    std::string HC = Sample.substr(0, 3) + "_" + Sample.substr(3, 4);

    TFile* Reg01 = TFile::Open((Sample + "Reg/PostPredRegVal0.1" + Sample + ".root").c_str()); // Load in histograms
    TFile* Reg05 = TFile::Open((Sample + "Reg/PostPredRegVal0.5" + Sample + ".root").c_str());
    TFile* Reg10 = TFile::Open((Sample + "Reg/PostPredRegVal1.0" + Sample + ".root").c_str()); // Load in histograms
    TFile* RegTrue01 = TFile::Open((Sample + "Reg/PostPredRegVal0.1" + Sample + "True.root").c_str()); // Load in histograms
    TFile* OscData = TFile::Open("EventRates/OscPMNSNoNC.root");
    TFile* RegNB = TFile::Open((Sample + "Reg/PostPredRegVal0.1" + Sample + "SB.root").c_str());

    auto HistoReg01 = Reg01->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str()); // Get titles of individual histograms
    auto HistoReg05 = Reg05->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str()); // Get titles of individual histograms
    auto HistoReg10 = Reg10->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str()); // Get titles of individual histograms
    auto HistoRegTrue01 = RegTrue01->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str()); // Get titles of individual histograms
    auto HistoData = OscData->Get<TH1D>(("hRecoNeutrinoEnergy" + HC).c_str());
    auto HistoRegNB = RegNB->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str());

    HistoReg01->SetLineColor(kCyan+3); // Cosmetics for oscillated histogram
    HistoReg01->SetLineWidth(2);
    HistoReg01->SetStats(0);

    HistoReg05->SetLineColor(kBlue); // Cosmetics for unoscillated histogram
    HistoReg05->SetLineWidth(2);
    HistoReg05->SetStats(0);

    HistoReg10->SetLineColor(kOrange+1); // Cosmetics for oscillated histogram
    HistoReg10->SetLineWidth(2);
    HistoReg10->SetStats(0);

    HistoRegTrue01->SetLineColor(kGreen+1); // Cosmetics for true histogram
    HistoRegTrue01->SetLineWidth(2);
    HistoRegTrue01->SetStats(0);

    HistoRegNB->SetLineColor(kRed); // Cosmetics for NB histogram
    HistoRegNB->SetLineWidth(2);
    HistoRegNB->SetStats(0);

    HistoData->SetLineColor(kBlack); // Cosmetics for data histogram
    HistoData->SetLineWidth(2);
    //HistoData->SetLineStyle(2);
    HistoData->SetStats(0);

    TCanvas *c = new TCanvas("c","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
    HistoRegNB->Draw("HIST");
    HistoRegNB->SetTitle(("Regularisation Strength Comparison - " + HC).c_str());
    //HistoReg05->Draw("HIST SAME");
    //HistoReg10->Draw("HIST SAME");
    //HistoRegTrue01->Draw("HIST SAME");
    HistoReg01->Draw("HIST SAME");
    HistoData->Draw("HIST SAME");

    auto leg = new TLegend(0.5, 0.5, 0.85, 0.85); // Creating the legend
    leg->SetMargin(0.4);
    
    leg->AddEntry(HistoRegNB,fmt::format("Reg 0.1 [0, 1] - Int = {:.2f}", HistoRegNB->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoReg01,fmt::format("Reg 0.1 [0, 2] - Int = {:.2f}", HistoReg01->Integral()).c_str(),"lpf");
    
    //leg->AddEntry(HistoReg05,fmt::format("Reg 0.5 - Int = {:.2f}", HistoReg05->Integral()).c_str(),"lpf");
    //leg->AddEntry(HistoReg10,fmt::format("Reg 1.0 - Int = {:.2f}", HistoReg10->Integral()).c_str(),"lpf");
    //leg->AddEntry(HistoRegTrue01,fmt::format("True 0.1 - Int = {:.2f}", HistoRegTrue01->Integral()).c_str(),"lpf");
    
    leg->AddEntry(HistoData,fmt::format("Data - Int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
    leg->Draw();

    c->SaveAs((Sample + "Reg/" + HC + "RegStrengthComp.png").c_str()); // Saving and creating title
}

void HPDMeanDist() {
    TFile* PostFile = TFile::Open("FHCnumuReg/RegValAll0.1_Process.root");
    auto HPDVal = PostFile->Get<TVectorT<double>>("Means_HPD");
    auto HPDErr = PostFile->Get<TVectorT<double>>("Errors_HPD");

    for (int i = 0; i < HPDVal->GetNrows(); i++) {
        std::cout << "HPD Mean for parameter " << i << " is: " << (*HPDVal)[i] << " with error: " << (*HPDErr)[i] << std::endl;
    }
}

void overlayHPD() { // Create a script to print two histograms on same canvas
    TFile* OscData = TFile::Open("EventRates/CCIndChanOsc.root");
    TFile* HPDER = TFile::Open("EventRates/FHCnumuIndChanHPD.root");

    auto HistoData = OscData->Get<TH1D>("FHC_numu_FHC_numubar_x_numubar"); // Get titles of individual histograms
    auto HistoHPD = HPDER->Get<TH1D>("FD_FHC_numu_FHC_numubar_x_numubar"); // Get titles of individual histograms

    HistoData->SetLineColor(kBlack); // Cosmetics for oscillated histogram
    HistoData->SetLineWidth(2);
    HistoData->SetStats(0);

    HistoHPD->SetLineColor(kBlue); // Cosmetics for unoscillated histogram
    HistoHPD->SetLineWidth(2);
    HistoHPD->SetStats(0);

    TCanvas *c = new TCanvas("c","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
    HistoHPD->Draw("HIST");
    HistoData->Draw("HIST SAME");

    auto leg = new TLegend(0.6, 0.7, 0.88, 0.88); // Creating the legend
    leg->AddEntry(HistoData,fmt::format("Data - Int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoHPD,fmt::format("HPD - Int = {:.2f}", HistoHPD->Integral()).c_str(),"l");
    leg->Draw();
    c->SaveAs("FHCnumuCompmubars.png"); // Saving and creating title
}

void RegDirComp(TString Sam) { // Create a script to print two histograms on same canvas
    std::string Sample = Sam.Data();
    std::string HC = Sample.substr(0, 3) + "_" + Sample.substr(3, 4);

    TFile* Reg = TFile::Open((Sample + "Reg/PostPredRegVal1.0" + Sample + "LB.root").c_str()); 
    TFile* RegL = TFile::Open(("PostPredHBinRSam.root")); // Load in histograms
    TFile* RegR = TFile::Open(("PostPredQBinRSam.root"));
    TFile* RegB = TFile::Open((Sample + "Reg/PostPredTrueCut20RegVal1" + Sample + ".root").c_str()); // Load in histograms
    TFile* OscData = TFile::Open("EventRates/OscPMNSNoNC.root");

    auto HistoReg = Reg->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str()); // Get titles of individual histograms
    auto HistoRegL = RegL->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str()); // Get titles of individual histograms
    auto HistoRegR = RegR->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str());
    auto HistoRegB = RegB->Get<TH1D>(("Predictive/FD_" + HC + "/FD_" + HC + "_mc_PostPred").c_str());
    auto HistoData = OscData->Get<TH1D>(("hRecoNeutrinoEnergy" + HC).c_str());
    
    HistoReg->SetLineColor(kOrange-3); // Cosmetics for oscillated histogram
    HistoReg->SetLineWidth(2);
    HistoReg->SetStats(0);
    
    HistoRegL->SetLineColor(kRed); // Cosmetics for oscillated histogram
    HistoRegL->SetLineWidth(2);
    HistoRegL->SetStats(0);

    HistoRegR->SetLineColor(kBlue); // Cosmetics for unoscillated histogram
    HistoRegR->SetLineWidth(2);
    HistoRegR->SetStats(0);

    HistoRegB->SetLineColor(kOrange-3); // Cosmetics for oscillated histogram
    HistoRegB->SetLineWidth(2);
    HistoRegB->SetStats(0);

    HistoData->SetLineColor(kBlack); // Cosmetics for data histogram
    HistoData->SetLineWidth(2);
    HistoData->SetLineStyle(2);
    HistoData->SetStats(0);

    TCanvas *c = new TCanvas("c","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
    HistoData->Draw("HIST");
    HistoData->SetTitle(("Parameter per Bin Comparison (R = 1.0) - " + HC).c_str());
    HistoRegL->Draw("HIST SAME");
    HistoRegR->Draw("HIST SAME");
    //HistoRegB->Draw("HIST SAME");
    HistoReg->Draw("HIST SAME");

    auto leg = new TLegend(0.5, 0.5, 0.85, 0.85); // Creating the legend
    leg->SetMargin(0.4);
    leg->AddEntry(HistoReg,fmt::format("1 Par/Bin - Int = {:.2f}", HistoReg->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoRegL,fmt::format("2 Par/Bin - Int = {:.2f}", HistoRegL->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoRegR,fmt::format("4 Par/Bin - Int = {:.2f}", HistoRegR->Integral()).c_str(),"lpf");
    //leg->AddEntry(HistoRegB,fmt::format("True Cut [0, 20] - Int = {:.2f}", HistoRegB->Integral()).c_str(),"lpf");
    leg->AddEntry(HistoData,fmt::format("Data - Int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
    leg->Draw();

    c->SaveAs((Sample + "Reg/" + HC + "ParPerBinComp.png").c_str()); // Saving and creating title
}

void CorrColPlot(TString file, int Chan){
    // 0 for numu->numu, 2 for numubar->numubar, 4 for numu->nue, 5 for numubar->nuebar
    TFile* CorrFile = TFile::Open(file);
    std::string FilePath = file.Data();
    FilePath.erase(FilePath.size() - 14);
    auto CorrPlot = CorrFile->Get<TH2D>("Correlation_plot");
    TCanvas *c = new TCanvas("hCorr","Correlation Plot");
    gStyle->SetPalette(104);
    int StartBin = (Chan * 40) + 1;
    int EndBin = (Chan + 1) * 40;
    std::string ChanName = CorrPlot->GetXaxis()->GetBinLabel(StartBin);
    ChanName.erase(ChanName.size() - 10);
    CorrPlot->GetXaxis()->SetRange(StartBin, EndBin);
    CorrPlot->GetYaxis()->SetRange(StartBin, EndBin);
    CorrPlot->SetTitle(("Correlation Plot - Channel " + ChanName).c_str());
    CorrPlot->Draw("COLZ");
    c->SaveAs((FilePath + "(" + ChanName + ")CorrCol.root").c_str());
    c->Close();
}