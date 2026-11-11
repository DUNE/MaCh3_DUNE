#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TImage.h"
#include "TCanvas.h"
#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <fmt/core.h>
#include <string>
#include <TH2D.h>
#include <iostream>
#include <vector>
#include <utility>
#include <map>
#include <limits>

// ############### POSTERIOR PREDICTIVE ANALYSIS ############### //

void PostPred(TString PPFile) {
    std::string FileName = PPFile.Data();
    std::string PlotFile = FileName.substr(0, FileName.size() - 5);
    std::string PDFName = PlotFile + "AllInfo.pdf";
    
    TFile* OscPMNS = TFile::Open("AllSampleFits/PostPredPMNS.root"); // Load in histograms
    TFile* OscMF = TFile::Open(PPFile);
    TFile* OscData = TFile::Open("EventRates/AllSampleTwentyOsc.root");
    TFile* OscFluc = TFile::Open("FlucData/FluctuatedData_5.root");

    TFile *OutFile = new TFile((PlotFile + "AllSampleInfo.root").c_str(), "RECREATE");
    std::vector<std::string> Samples = {"FHC_numu", "FHC_nue", "RHC_numu", "RHC_nue"};
    bool firstPage = true;
    
    for (const auto& Sample : Samples) { 

        auto HistoPMNS = OscPMNS->Get<TH1D>(("Predictive/" + Sample + "/" + Sample + "_mc_PostPred").c_str()); // PostPredPMNSv2
        auto HistoMF = OscMF->Get<TH1D>(("Predictive/FD_" + Sample + "/FD_" + Sample + "_mc_PostPred").c_str()); // PostPredMFv2
        auto HistoDim0 = OscMF->Get<TH2D>(("Predictive/FD_" + Sample + "/FD_" + Sample + "_mc_dim0").c_str());
        auto HistoData = OscData->Get<TH1D>(("hFD_" + Sample).c_str());
        auto HistoFluc = OscFluc->Get<TH1D>(("hFD_" + Sample).c_str());

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

        HistoFluc->SetLineColor(kBlack); // Cosmetics for data histogram
        HistoFluc->SetLineWidth(2);
        HistoFluc->SetLineStyle(2);
        HistoFluc->SetStats(0);

        TCanvas *c = new TCanvas(("PostPred" + Sample).c_str(),"Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 

        HistoDim0->Draw("COLZ");
        //HistoPMNS->Draw("HIST E0 SAME");
        HistoMF->Draw("HIST E0 SAME");

        if (PlotFile.find("Fluc") != std::string::npos){
            HistoFluc->Draw("hist same");
        } 
        else {
            HistoData->Draw("HIST SAME");
        }

        HistoDim0->SetTitle(("Posterior Predictive - " + Sample).c_str());
        HistoDim0->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
        HistoDim0->GetYaxis()->SetTitle("Number of Events");
        HistoDim0->GetZaxis()->SetTitle("Frequency of reweighted toy event counts");

        c->SetLogz();
        c->SetRightMargin(0.15);
        gStyle->SetPalette(kCividis);

        HistoMF->GetXaxis()->SetRangeUser(0, 10);
        HistoPMNS->GetXaxis()->SetRangeUser(0, 10);
        HistoDim0->GetXaxis()->SetRangeUser(0, 10);
        HistoData->GetXaxis()->SetRangeUser(0, 10);
        HistoFluc->GetXaxis()->SetRangeUser(0, 10);

        auto leg = new TLegend(0.35, 0.5, 0.84, 0.84); // Creating the legend
        leg->SetMargin(0.4);
        leg->AddEntry(HistoMF,fmt::format("Model-light - #int = {:.2f}", HistoMF->Integral()).c_str(),"lpf");
        //leg->AddEntry(HistoPMNS,fmt::format("PMNS - #int = {:.2f}", HistoPMNS->Integral()).c_str(),"lpf");
        if (PlotFile.find("Fluc") != std::string::npos){
            leg->AddEntry(HistoFluc,fmt::format("Fluc Data - #int = {:.2f}", HistoFluc->Integral()).c_str(), "lpf");
        }
        else {
            leg->AddEntry(HistoData,fmt::format("Data - #int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
        }
        leg->Draw();

        OutFile->cd();
        if (firstPage){
            c->Print((PDFName + "(").c_str());
            firstPage = false;
        }
        else{
            c->Print(PDFName.c_str());
        }
        c->Write();

        TCanvas *g = nullptr;
        OscMF->GetObject(("Predictive/FD_" + Sample + "/FD_" + Sample + "_sum_canv").c_str(), g);
        g->Draw();
        g->SetName(("SumCanv" + Sample).c_str());
        OutFile->cd();
        g->Write();
        g->Print(PDFName.c_str());

        TCanvas *h = nullptr;
        OscMF->GetObject(("Predictive/FD_" + Sample + "/FD_" + Sample + "_predfluc_draw_canv").c_str(), h);
        h->Draw();
        h->SetName(("PredFlucDraw" + Sample).c_str());
        OutFile->cd();
        h->Write();
        h->Print(PDFName.c_str());

        TCanvas *j = nullptr;
        OscMF->GetObject(("Predictive/FD_" + Sample + "/FD_" + Sample + "_drawfluc_draw_canv").c_str(), j);
        j->Draw();
        j->SetName(("DrawFlucDraw" + Sample).c_str());
        OutFile->cd();
        j->Write();
        j->Print(PDFName.c_str());

    }
    OutFile->Close();
    TCanvas *closer = new TCanvas();
    closer->Print((PDFName + ")").c_str());
}

void PMNSParVal(){
    TFile* DataOsc = TFile::Open("EventRates/AllSampleTwentyTrueIndChanOsc.root");
    TFile* DataUnosc = TFile::Open("EventRates/AllSampleTwentyTrueIndChanUnosc.root");
    TFile* ColMap = TFile::Open("AllSampleTwentyTestOscProbCol.root");

    std::vector<std::string> Samples = {"FHC_numu", "FHC_nue", "RHC_numu", "RHC_nue"};
    std::vector<std::string> Channel = {"FD_FHC_numu_FHC_numu_x_numu", "FD_FHC_nue_FHC_numu_x_nue", "FD_RHC_numu_RHC_numubar_x_numubar", "FD_RHC_nue_RHC_numubar_x_nuebar"};
    std::vector<std::string> OscMode = {"numu_to_numu", "numu_to_nue", "numubar_to_numubar", "numubar_to_nuebar"};

    TFile* outFile = new TFile("RatioParam.root", "RECREATE");

    for (int i = 0; i < 4; ++i) {     
        auto DataO = DataOsc->Get<TH1D>(Channel[i].c_str());
        auto DataU = DataUnosc->Get<TH1D>(Channel[i].c_str());

        TH1D* DataR = (TH1D*)DataO->Clone(("DataRat_" + Samples[i]).c_str());
        for (int i = 1; i <= DataO->GetNbinsX(); i++){
            double Oval = DataO->GetBinContent(i);
            double Uval = DataU->GetBinContent(i);
            if (Uval == 0) DataR->SetBinContent(i, 0);
            else DataR->SetBinContent(i, Oval / Uval);
        }

        TCanvas *c = new TCanvas(("PostPred" + Samples[i]).c_str(),"Overlay");
        ColMap->GetObject(OscMode[i].c_str(), c);
        c->Draw();
        DataR->SetLineColor(kBlack);
        DataR->Draw("hist same");
        TLegend* leg = new TLegend(0.4, 0.75, 0.6, 0.85);
        leg->AddEntry(DataR, "PMNS Probability", "l");
        leg->Draw();
        c->Write();
    }

    outFile->Close();
}

// ############## PARAMETER POSTERIOR DISTRIBUTION ############## //

const int NParams = 960; // Number of energy normalisation parameters

const int NEBins = 80; // Number of energy bins per channel
const double EMin = 0.0;
const double EMax = 20.0;
const double BinWidth = (EMax - EMin) / NEBins;

const int NProbBins = 1000; // Number of probability bins in the colour plots
const double ProbMin = 0.0;
const double ProbMax = 2.0;
const double ProbBinWidth = (ProbMax - ProbMin)/ NProbBins;

std::vector<double> NormParams(NParams, 0.0); // Creating empty vector to store all energy normalisation parameters
std::vector<TH2D*> EventHistos; // Creating empty vector to store histograms
std::vector<double> BinLowEdges;
std::vector<double> ProbBinLowEdges;
std::vector<double> BinCentres;
std::vector<std::string> OscChans = {"numu_to_numu", "nue_to_nue", "numubar_to_numubar", "nuebar_to_nuebar", "numu_to_nue", "numubar_to_nuebar", "nue_to_numu", "nuebar_to_numubar", "numu_to_nutau", "nue_to_nutau", "numubar_to_nutaubar", "nuebar_to_nutaubar"};

void SetupBinning(){
    for (int a = 0; a <= NEBins; a++){
        double BinLow = EMin + a * BinWidth;
        BinLowEdges.push_back(BinLow);
    }
    for (int b = 0; b < NEBins; b++){
        double BinMid = (BinLowEdges[b] + BinLowEdges[b+1]) / 2;
        BinCentres.push_back(BinMid);
    }
    for (int c = 0; c <= NProbBins; c++){
        double ProbBinLow = ProbMin + c * ProbBinWidth;
        ProbBinLowEdges.push_back(ProbBinLow);
    }
}

void NormaliseAndPlot(TH2D* Histo, TString HistoName, int PlotEntries, TString File_To_Plot) // Create a script to column normalise and then plot a given histogram with a given name
{

    for (int x = 1; x <= Histo->GetNbinsX(); x++) { // For all energy bins, avoiding underflow bin
        bool FixedCol = true; // Assume a column is for a fixed parameter to start
        for (int y = 2; y <= Histo->GetNbinsY(); y++) { // For every bin except the first (if fixed, all events will be in bin 1)
            if (Histo->GetBinContent(x,y) != 0.0) { // If bin content is not zero
                FixedCol = false; // Not a fixed parameter - only need one non-zero bin at not zero to prove this
                break;
            }
        }
        if (FixedCol) { // If not fixed then skip this, if fixed then;
            for (int y = 1; y <= Histo->GetNbinsY(); y++) {
                Histo->SetBinContent(x, y, 0.0); // Set all bin contents to be NaN, will all show up dark blue
            }
            continue; 
        }
        double ColSum = 0;
        for (int y = 1; y <= Histo->GetNbinsY(); y++) { // For all oscillation probability bins, avoiding underflow bin
            ColSum += Histo->GetBinContent(x,y); // Running total of the number of events in the column
        }
        if (ColSum == 0) continue;
        for (int y = 1; y <= Histo->GetNbinsY(); y++) { // For all oscillation probability bins, avoiding underflow bin
            double Content = Histo->GetBinContent(x,y); // Get the content of each bin
            Histo->SetBinContent(x, y, Content/ColSum); // Normalise each bin by the total number of events in the column
        }
    }

    bool isMuChan = false;
    bool isEChan = false;
    if (HistoName.EndsWith("mu") || HistoName.EndsWith("mubar")) isMuChan = true;
    if (HistoName.EndsWith("e") || HistoName.EndsWith("ebar")) isEChan = true;

    std::string HistNam = HistoName.Data();
    std::string FileName = File_To_Plot.Data();
    std::string PlotFile = FileName.substr(0, FileName.size() - 5);
    TCanvas *Canvas = new TCanvas(HistNam.c_str(),"Overlay"); // Create a canvas to plot
    Canvas->Draw(); // Draw the canvas
    Canvas->cd();
    if(isMuChan) Canvas->SetLogy();
    Histo->Draw("COLZ"); // Draw colour scale histogram
    Histo->GetXaxis()->SetTitle("True Neutrino Energy (GeV)"); // Setting axes properties
    Histo->GetXaxis()->SetRangeUser(0, 10);
    Histo->GetYaxis()->SetTitle("Oscillation Probability");
    if(isEChan) Histo->GetYaxis()->SetRangeUser(0, 0.2);
    Histo->GetZaxis()->SetTitle(("Number of Events (Total = " + std::to_string(PlotEntries) + ")").c_str());
    Histo->GetZaxis()->SetTitleOffset(1.5);
    Histo->SetStats(0); // Cosmetics of plot
    Canvas->SetRightMargin(0.15);
    Canvas->Update();
    gDirectory->cd();
    Canvas->Write();
}

// ############## LOADING BIN SCHEME, DATA, AND PLOTTING ############## //

void ChanColPlot(TString File_To_Plot)//, int Chan_To_Plot) // Create a script to produce colour plots for each oscillation channel from a given root file
{
    BinLowEdges.clear();
    ProbBinLowEdges.clear();
    BinCentres.clear();
    EventHistos.clear();

    std::string FileName = File_To_Plot.Data();
    std::string PlotFile = FileName.substr(0, FileName.size() - 5);

    SetupBinning();

    auto RootFile = TFile::Open(File_To_Plot); // Open the specified file

    if (RootFile->IsZombie() || !RootFile) { // Check if file exists
        throw std::runtime_error("File does not exist");
    }

    auto Posteriors = RootFile->Get<TTree>("posteriors"); // Read in all posteriors from file
    int numTotPars = Posteriors->GetListOfLeaves()->GetEntries(); 
    int parsRunOver = numTotPars - 7;

    int parIdx = 0;
    //for (int i = 0; i < 960; i++) {
    for (int i = (parsRunOver - 960); i < parsRunOver; i++) { // For all energy normalisation parameters
        std::string ChainVarName = "param_" + std::to_string(i); // Create relevent variable name
        Posteriors->SetBranchAddress(ChainVarName.c_str(), &NormParams[parIdx]); // Reassign address to correct name and store values in given vector
        parIdx++;
    }

    for (const auto& channel : OscChans) { // For each channel in our pairing of oscillation channels
        EventHistos.push_back(new TH2D(channel.c_str(), ("Oscillation Probability for " + channel).c_str(), NEBins, BinLowEdges.data(), NProbBins, ProbBinLowEdges.data())); // Create histogram for each channel in each sample
    }

    Long64_t MCEntries = Posteriors->GetEntries(); // Get number of entries (steps) in given file

    const int BurnInCut = 0.25 * MCEntries; // Number of steps to cut out for burn in
    const int PlotEntries = MCEntries - BurnInCut; // Number of steps to plot after burn in

    for (Long64_t i = BurnInCut; i < MCEntries; i++) { // Start at some burn-in cut, run for rest of entries
        Posteriors->GetEntry(i); // Get the relevant information
        for (int p = 0; p < 960; p++) { // For each parameter
            int k = p % NEBins; // Find the value of param index modulo 40 (each value corresponds to a different energy bin (40))
            int g = p / NEBins; // Find the integer division of param index / 48 (each value corresponds to a different histogram (48))
            EventHistos[g]->Fill(BinCentres[k], NormParams[p]); // Fill each histogram with the relevant parameter from the desired energy bin
        }
    }

    std::vector<TH2D*> MajorChan = {EventHistos[0], EventHistos[2], EventHistos[4], EventHistos[5]};

    TFile* outFile = new TFile((PlotFile + "OscProbCol.root").c_str(), "RECREATE");

    for (const auto& Plot : MajorChan){
        NormaliseAndPlot(Plot, Plot->GetName(), PlotEntries, File_To_Plot);
    }

    outFile->Close();

}


void OverlapPostPred(TString PPFile1, TString PPFile2) {

    std::string PDFName = "PostPredComp.pdf";
    
    TFile* OscPMNS = TFile::Open("AllSampleFits/PostPredPMNS.root"); // Load in histograms
    TFile* OscMF = TFile::Open(PPFile1);
    TFile* OscHenry = TFile::Open(PPFile2);
    TFile* OscData = TFile::Open("EventRates/AllSampleTwentyOsc.root");

    TFile *OutFile = new TFile("HenryCompAllSampleInfo.root", "RECREATE");
    std::vector<std::string> Samples = {"FHC_numu", "FHC_nue", "RHC_numu", "RHC_nue"};
    bool firstPage = true;
    
    for (const auto& Sample : Samples) { 

        auto HistoPMNS = OscPMNS->Get<TH1D>(("Predictive/" + Sample + "/" + Sample + "_mc_PostPred").c_str()); // PostPredPMNSv2
        auto HistoMF = OscMF->Get<TH1D>(("Predictive/FD_" + Sample + "/FD_" + Sample + "_mc_PostPred").c_str()); // PostPredMFv2
        auto HistoData = OscData->Get<TH1D>(("hFD_" + Sample).c_str());
        auto HistoFluc = OscHenry->Get<TH1D>(("Predictive/FD_" + Sample + "/FD_" + Sample + "_mc_PostPred").c_str());

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

        HistoFluc->SetLineColor(kGreen); // Cosmetics for data histogram
        HistoFluc->SetLineWidth(2);
        HistoFluc->SetStats(0);

        TCanvas *c = new TCanvas(("PostPred" + Sample).c_str(),"Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 

        HistoMF->Draw("HIST E0");
        HistoFluc->Draw("hist same");
        HistoData->Draw("HIST SAME");
        HistoPMNS->Draw("hist same");

        HistoMF->SetTitle(("Posterior Predictive - " + Sample).c_str());
        HistoMF->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
        HistoMF->GetYaxis()->SetTitle("Number of Events");
        HistoMF->GetZaxis()->SetTitle("Frequency of reweighted toy event counts");

        c->SetLogz();
        c->SetRightMargin(0.15);

        HistoMF->GetXaxis()->SetRangeUser(0, 10);
        HistoPMNS->GetXaxis()->SetRangeUser(0, 10);
        HistoData->GetXaxis()->SetRangeUser(0, 10);
        HistoFluc->GetXaxis()->SetRangeUser(0, 10);

        auto leg = new TLegend(0.4, 0.6, 0.84, 0.84); // Creating the legend
        leg->SetMargin(0.4);
        leg->AddEntry(HistoMF,fmt::format("Max 2 - #int = {:.2f}", HistoMF->Integral()).c_str(),"lpf");
        leg->AddEntry(HistoPMNS,fmt::format("PMNS - #int = {:.2f}", HistoPMNS->Integral()).c_str(),"lpf");
        leg->AddEntry(HistoFluc,fmt::format("Max 1 - #int = {:.2f}", HistoFluc->Integral()).c_str(), "lpf");
        leg->AddEntry(HistoData,fmt::format("Data - #int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
        leg->Draw();

        OutFile->cd();
        if (firstPage){
            c->Print("NewPDF.pdf(");
            firstPage = false;
        }
        else{
            c->Print("NewPDF.pdf");
        }
        c->Write();

    }
    OutFile->Close();
    TCanvas *closer = new TCanvas();
    closer->Print("NewPDF.pdf)");
}