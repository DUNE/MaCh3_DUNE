#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TImage.h"
#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <fmt/core.h>
#include <unordered_map>
#include <string>
#include <iostream>

void PoissonFluc(){
    TFile* LoadData = TFile::Open("EventRates/AllSampleTwentyOsc.root");
    std::vector<std::string> Histos = {"hFD_FHC_numu", "hFD_FHC_nue", "hFD_RHC_numu", "hFD_RHC_nue"};
    TRandom3 rng(0);
    for (int flucNum = 1; flucNum <= 10; flucNum++){
        std::string outName = "FlucData/FluctuatedData_" + std::to_string(flucNum) + ".root";
        TFile* flucData = TFile::Open(outName.c_str(), "RECREATE");
        for (const auto& samp : Histos){
            auto dataHist = LoadData->Get<TH1D>(samp.c_str());
            auto newHist = dynamic_cast<TH1*>(dataHist->Clone());
            newHist->Reset();
            for (int i = 1; i <= dataHist->GetNbinsX(); i++){
                double events = dataHist->GetBinContent(i);
                double flucEvents = rng.Poisson(events);
                newHist->SetBinContent(i, flucEvents);
            }
            newHist->Write();
        }
        flucData->Close();
    }
    LoadData->Close();
}

void overlay() { // Create a script to print two histograms on same canvas
    TFile* OscMF = TFile::Open("EventRates/AllSampleTwentyOsc.root"); // Load in histograms
    TFile* OscPMNS = TFile::Open("RecoTrueER.root");
    TFile* Extra = TFile::Open("RecoERFHCnue.root");

    TFile *OutFile = new TFile("TruevPMNSER_AllSamples.root", "RECREATE");
    TIter next(OscMF->GetListOfKeys()); // Producing list of keys (histograms) from file
    TKey* key; 
    while ((key = (TKey*)next())) { 
        std::string KeyName = key->GetName();
        auto HistoOsc = OscMF->Get<TH1D>(KeyName.c_str()); // Get titles of individual histograms
        auto HistoUnosc = OscPMNS->Get<TH1D>(KeyName.c_str()); // Get titles of individual histograms
        auto HistoExtra = Extra->Get<TH1D>(KeyName.c_str());

        HistoOsc->SetLineColor(kBlack); // Cosmetics for oscillated histogram
        HistoOsc->SetLineWidth(2);
        HistoOsc->SetStats(0);

        HistoUnosc->SetLineColor(kMagenta); // Cosmetics for unoscillated histogram
        HistoUnosc->SetLineWidth(2);
        HistoUnosc->SetStats(0);

        HistoExtra->SetLineColor(kOrange); // Cosmetics for unoscillated histogram
        HistoExtra->SetLineWidth(2);
        HistoExtra->SetStats(0);

        TCanvas *c = new TCanvas(KeyName.c_str(),"Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
        HistoOsc->Draw("HIST");
        HistoOsc->SetTitle(("Event Rate Comparison - " + KeyName).c_str());
        HistoOsc->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
        HistoOsc->GetYaxis()->SetTitle("Number of Events");
        HistoUnosc->Draw("HIST SAME");
        HistoExtra->Draw("HIST SAME");

        auto leg = new TLegend(0.6, 0.6, 0.88, 0.88); // Creating the legend
        leg->AddEntry(HistoUnosc,fmt::format("True Params - #int = {:.2f}", HistoUnosc->Integral()).c_str(),"lpf");
        leg->AddEntry(HistoOsc,fmt::format("Data - #int = {:.2f}", HistoOsc->Integral(0, 40)).c_str(),"l");
        leg->AddEntry(HistoExtra,fmt::format("Reco Params - #int = {:.2f}", HistoOsc->Integral()).c_str(),"l");
        leg->Draw();
        OutFile->cd();
        c->Write();
        //c->SaveAs(("TruevPMNSER(" + Sample + ").root").c_str()); // Saving and creating title
    }
    OutFile->Close();
}

void PostPred(TString PPFile) { // Create a script to print two histograms on same canvas
    std::string FileName = PPFile.Data();
    std::string PlotFile = FileName.substr(0, FileName.size() - 5);
    std::string PDFName = PlotFile + "AllInfo.pdf";
    
    TFile* OscPMNS = TFile::Open("AllSampleFits/PostPredPMNS.root"); // Load in histograms
    TFile* OscMF = TFile::Open(PPFile);
    TFile* OscData = TFile::Open("EventRates/AllSampleTwentyOsc.root");
    TFile* OscFluc = TFile::Open("FlucData/FluctuatedData_1.root");

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
        HistoPMNS->Draw("HIST E0 SAME");
        HistoMF->Draw("HIST E0 SAME");
        HistoData->Draw("HIST SAME");
        HistoFluc->Draw("hist same");
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
        //HistoDim0->GetYaxis()->SetRangeUser(0, 40);

        auto leg = new TLegend(0.35, 0.5, 0.84, 0.84); // Creating the legend
        leg->SetMargin(0.4);
        leg->AddEntry(HistoMF,fmt::format("Model-light - #int = {:.2f}", HistoMF->Integral()).c_str(),"lpf");
        leg->AddEntry(HistoPMNS,fmt::format("PMNS - #int = {:.2f}", HistoPMNS->Integral()).c_str(),"lpf");
        leg->AddEntry(HistoData,fmt::format("Data - #int = {:.2f}", HistoData->Integral()).c_str(),"lpf");
        leg->AddEntry(HistoFluc,fmt::format("Fluc Data - #int = {:.2f}", HistoFluc->Integral()).c_str(), "lpf");
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

// Tidy up this macro tomorrow!

void ToysPostComp(){
    auto ParamTree = TFile::Open("FHCnueReg/PostPredRegVal0.1FHCnue.root")->Get<TTree>("ToySummary");
    auto Posterior = TFile::Open("FHCnueReg/RegVal0.1FHCnue_Process.root")->Get<TH1D>(("Post_1d_hists/numu_nue_2.50_2.75"));
    TFile* PostFile = TFile::Open("FHCnueReg/RegVal0.1FHCnue_Process.root");
    auto HPDVal = PostFile->Get<TVectorT<double>>("Means_HPD");
    auto HPD = (*HPDVal)[170];
    std::cout << "HPD Mean for Param 170 is " << HPD << std::endl;
    auto PDFVal = PostFile->Get<TVectorT<double>>("PDF_Means");
    auto PDF = (*PDFVal)[170];
    std::cout << "PDF Mean for Param 170 is " << PDF << std::endl;
    auto GaussVal = PostFile->Get<TVectorT<double>>("Gauss_Means");
    auto Gauss = (*GaussVal)[170];
    std::cout << "Gauss Mean for Param 170 is " << Gauss << std::endl;
    double value = 0.0;
    ParamTree->SetBranchAddress("numu_nue_2.50_2.75", &value);
    TH1D *h = new TH1D("hToy", "Toy Posterior Distribution", Posterior->GetXaxis()->GetNbins(), Posterior->GetXaxis()->GetXmin(), Posterior->GetXaxis()->GetXmax());
    Long64_t nEntries = ParamTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        ParamTree->GetEntry(i);
        h->Fill(value);
    }
    TCanvas *c = new TCanvas("ToysAgainstPosterior","Toy Posterior Distribution");
    h->Scale(1.0 / h->Integral()); // Normalize the histogram
    Posterior->Scale(1.0 / Posterior->Integral()); // Normalize the histogram
    int maxBin = h->GetMaximumBin();
    h->Draw("HIST");
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->SetStats(0);
    std::cout << "Toy Mean for Param 192 is " << h->GetMean() << std::endl;
    h->GetXaxis()->SetRangeUser(0, 0.2);
    Posterior->Draw("HIST SAME");
    //h->GetYaxis()->SetRangeUser(0, Posterior->GetMaximum() * 1.1);
    Posterior->SetLineColor(kRed);
    Posterior->SetLineWidth(2);
    auto leg = new TLegend(0.5, 0.5, 0.88, 0.88);
    leg->AddEntry(h,fmt::format("Toy Distribution - Mean = {:.4f}", h->GetMean()).c_str(),"lpf");
    leg->AddEntry(Posterior, fmt::format("Parameter Posterior - HPD = {:.4f}", HPD).c_str(),"l");
    leg->AddEntry(Posterior, fmt::format("Parameter Posterior - PDF = {:.4f}", PDF).c_str(),"l");
    leg->AddEntry(Posterior, fmt::format("Parameter Posterior - Gauss = {:.4f}", Gauss).c_str(),"l");
    leg->Draw();
    c->SaveAs(("ToyPostFHCnue(numu_nue_2.50_2.75).root"));
    c->SaveAs(("ToyPostFHCnue(numu_nue_2.50_2.75).png"));
}

void PlotSigmaVar(TString inputFile, TString paramName, TString sampleName){
    std::string fileName = inputFile.Data();
    std::string param = paramName.Data();
    std::string sample = sampleName.Data();
    TFile* sigmaFile = TFile::Open(inputFile);
    TFile *out = new TFile(("SigmaVar_FHC" + sample + "(" + param + ").root").c_str(), "RECREATE");
    TIter nextKey(sigmaFile->GetListOfKeys());
    TKey* key;
    bool firstPage = true;
    while ((key = (TKey*)nextKey())){
        TString dirName = key->GetName();
        if (!dirName.BeginsWith(paramName)) continue;
        TDirectoryFile *chanDir = sigmaFile->Get<TDirectoryFile>(dirName);
        if (!chanDir) continue;
        TIter nextSubKey(chanDir->GetListOfKeys());
        TKey* subKey;
        while ((subKey = (TKey*)nextSubKey())){
            TString sampleFolder = subKey->GetName();
            if (!sampleFolder.EndsWith(sampleName)) continue;
            TDirectoryFile *sampleDir = chanDir->Get<TDirectoryFile>(sampleFolder);
            if (!sampleDir) continue;
            std::vector<int> lineCols = {kRed, kGreen+1, kBlack, kBlue+1, kOrange+1};
            std::vector<int> lineStyles = {kDotted, kDashed, kSolid, kDashDotted, kDashDotted};
            std::vector<std::string> lineLabels = {"-3#sigma", "-1#sigma", "Nominal", "+1#sigma", "+3#sigma"};
            TCanvas *c = new TCanvas(dirName);
            auto leg = new TLegend(0.6, 0.6, 0.88, 0.88);
            double ymax = 0.0;
            for (int i = 0; i < 5; i++) {
                TH1D *h = sampleDir->Get<TH1D>(Form("Variation_%d", i));
                ymax = std::max(ymax, h->GetMaximum());
            }
            for (int i = 0; i < 5; i++){
                auto sigmaHist = sampleDir->Get<TH1D>(("Variation_" + std::to_string(i)).c_str());
                sigmaHist->SetLineColor(lineCols[i]);
                sigmaHist->SetLineStyle(lineStyles[i]);
                sigmaHist->SetLineWidth(2);
                sigmaHist->SetStats(0);
                leg->AddEntry(sigmaHist, lineLabels[i].c_str(), "l");
                if (i == 0){
                    sigmaHist->Draw("HIST");
                    TString histName = Form("SigmaVar_%s(%s)", sampleFolder.Data(), dirName.Data());
                    sigmaHist->SetTitle(histName);
                    sigmaHist->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
                    sigmaHist->GetYaxis()->SetTitle("Number of Events");
                    sigmaHist->SetMaximum(ymax * 1.1);
                } else {
                    sigmaHist->Draw("HIST SAME");
                }
            }
            leg->Draw();
            out->cd();
            c->Write();

            if (firstPage) {
                c->SaveAs(("SigmaVar_FHC" + sample + "(" + param + ").pdf(").c_str());   // open PDF
                firstPage = false;
            } else {
                c->SaveAs(("SigmaVar_FHC" + sample + "(" + param + ").pdf").c_str());    // append page
            }
        }
    }
    TCanvas *closer = new TCanvas();
    closer->SaveAs(("SigmaVar_FHC" + sample + "(" + param + ").pdf)").c_str());

    out->Close();
    sigmaFile->Close();
}

void PlotCorr(){
    TFile* CorrFile = TFile::Open("FHCnueAllChanFree.root");
    auto postTree = CorrFile->Get<TTree>("posteriors");
    TFile *outFile = new TFile("PosteriorFHCnueAllChanFree.root", "RECREATE");
    std::vector<int> idx;
    bool firstPage = true;
    for (int k = 0; k < 960; k++) {
        if (k % 80 == 25) idx.push_back(k);
    }
    for (int k = 0; k < 960; k++){
        if (k % 80 == 25){
            std::string Params = "param_" + std::to_string(k);
            TCanvas *c = new TCanvas(Params.c_str(), Params.c_str(), 800, 600);
            postTree->Draw(Params.c_str(), "step>5000000", "l");
            outFile->cd();
            c->Write();
            if (firstPage){
                c->SaveAs("PosteriorPlots.pdf(");
                firstPage = false;
            }
            else{
                c->SaveAs("PosteriorPlots.pdf");
            }
            delete c;
        }
    }
    for (int a = 0; a < idx.size(); a++){
        for (int b = a + 1; b < idx.size(); b++){
            int j = idx[a];
            int k = idx[b];
            std::string ParamPair = "param_" + std::to_string(j) + ":param_" + std::to_string(k);
            TCanvas *c = new TCanvas(ParamPair.c_str(), ParamPair.c_str(), 800, 600);
            postTree->Draw(ParamPair.c_str(), "step>5000000", "colz");
            outFile->cd();
            c->Write();
            c->SaveAs("PosteriorPlots.pdf");
            delete c;
        }
    }
    TCanvas *closer = new TCanvas();
    closer->SaveAs("PosteriorPlots.pdf)");
    delete closer;
    outFile->Close();
    CorrFile->Close();
}

void drawCanvas2x2() {

    TFile *f = TFile::Open("PostPredAllSampleTwentyHalfAllSampleInfo.root", "READ");

    // Retrieve canvases by name
    TCanvas *c1 = nullptr;
    TCanvas *c2 = nullptr;
    TCanvas *c3 = nullptr;
    TCanvas *c4 = nullptr;

    f->GetObject("PostPredFHC_numu", c1);
    f->GetObject("PostPredFHC_nue", c2);
    f->GetObject("PostPredRHC_numu", c3);
    f->GetObject("PostPredRHC_nue", c4);

    // Null checks
    if (!c1) std::cerr << "Missing canvas: can1" << std::endl;
    if (!c2) std::cerr << "Missing canvas: can2" << std::endl;
    if (!c3) std::cerr << "Missing canvas: can3" << std::endl;
    if (!c4) std::cerr << "Missing canvas: can4" << std::endl;

    if (!c1 || !c2 || !c3 || !c4) {
        std::cerr << "ERROR: One or more canvases missing. Aborting." << std::endl;
        return;
    }

    const char* f1 = "AllSampleTwentyTest(numu_to_numu)OscProbCol.png";
    const char* f2 = "AllSampleTwentyTest(numu_to_nue)OscProbCol.png";
    const char* f3 = "AllSampleTwentyTest(numubar_to_numubar)OscProbCol.png";
    const char* f4 = "AllSampleTwentyTest(numubar_to_nuebar)OscProbCol.png";
    
    TFile *outFile = new TFile("4PlotParamPost.root", "RECREATE");

    // Create a new 2x2 canvas
    TCanvas *big = new TCanvas("big", "2x2 Canvas Grid", 1200, 1000);
    big->Divide(2, 2);

    // --- Pad 1 ---
    big->cd(1);
    c1->DrawClonePad();
    // TImage* img1 = TImage::Open(f1);
    // if (img1) img1->Draw("xxx");
    TLatex *lab1 = new TLatex(0.05, 0.95, "(a)");
    lab1->SetNDC();
    lab1->SetTextSize(0.06);
    lab1->Draw();

    // --- Pad 2 ---
    big->cd(2);
    c2->DrawClonePad();
    //TImage* img2 = TImage::Open(f2);
    //if (img2) img2->Draw("xxx");
    TLatex *lab2 = new TLatex(0.05, 0.95, "(b)");
    lab2->SetNDC();
    lab2->SetTextSize(0.06);
    lab2->Draw();

    // --- Pad 3 ---
    big->cd(3);
    c3->DrawClonePad();
    //TImage* img3 = TImage::Open(f3);
    //if (img3) img3->Draw("xxx");
    TLatex *lab3 = new TLatex(0.05, 0.95, "(c)");
    lab3->SetNDC();
    lab3->SetTextSize(0.06);
    lab3->Draw();

    // --- Pad 4 ---
    big->cd(4);
    c4->DrawClonePad();
    //TImage* img4 = TImage::Open(f4);
    //if (img4) img4->Draw("xxx");
    TLatex *lab4 = new TLatex(0.05, 0.95, "(d)");
    lab4->SetNDC();
    lab4->SetTextSize(0.06);
    lab4->Draw();

    big->Update();
    big->Write();
    big->SaveAs("4PlotParamPost.png");

    outFile->Close();
}

void RegCompPost(){
    TFile* BestReg = TFile::Open("FHCnumuReg/PostPredRegVal0.1FHCnumu.root"); // Load in histograms
    TFile* StrongestReg = TFile::Open("FHCnumuReg/PostPredRegVal20FHCnumu.root");
    TFile* NoReg = TFile::Open("FHCnumuOnly/PostPredFHCnumu.root");
    TFile* OscData = TFile::Open("EventRates/AllSampleTwentyOsc.root");

    auto BestHist = BestReg->Get<TH1D>("Predictive/FD_FHC_numu/FD_FHC_numu_mc_PostPred");
    auto StrongHist = StrongestReg->Get<TH1D>("Predictive/FD_FHC_numu/FD_FHC_numu_mc_PostPred");
    auto NoneHist = NoReg->Get<TH1D>("Predictive/FHC_numu/FHC_numu_mc_PostPred");
    auto HistoData = OscData->Get<TH1D>("hFD_FHC_numu");

    BestHist->SetLineColor(kBlue); // Cosmetics for PMNS model histogram
    BestHist->SetLineWidth(2);
    BestHist->SetStats(0);

    StrongHist->SetLineColor(kOrange); // Cosmetics for Model Free model histogram
    StrongHist->SetLineWidth(2);
    StrongHist->SetStats(0);

    NoneHist->SetLineColor(kViolet); // Cosmetics for Model Free model histogram
    NoneHist->SetLineWidth(2);
    NoneHist->SetStats(0);

    HistoData->SetLineColor(kBlack); // Cosmetics for data histogram
    HistoData->SetLineWidth(2);
    HistoData->SetLineStyle(2);
    HistoData->SetStats(0);

    TCanvas *c = new TCanvas("RegStr","Overlay"); // Creating canvas and plotting, may want to change order of Osc/Unosc based on visuals 
    BestHist->Draw("hist");
    StrongHist->Draw("hist same");
    NoneHist->Draw("hist same");
    HistoData->Draw("hist same");
    BestHist->SetTitle("Regularisation Comparison - FHC_numu");
    BestHist->GetXaxis()->SetTitle("Reconstructed Neutrino Energy (GeV)");
    BestHist->GetYaxis()->SetTitle("Number of Events");
    BestHist->GetXaxis()->SetRangeUser(0, 10);
    StrongHist->GetXaxis()->SetRangeUser(0, 10);
    NoneHist->GetXaxis()->SetRangeUser(0, 10);
    HistoData->GetXaxis()->SetRangeUser(0, 10);

    auto leg = new TLegend(0.45, 0.45, 0.85, 0.85); // Creating the legend
    leg->AddEntry(BestHist,"Reg Strength 0.1","l");
    leg->AddEntry(StrongHist,"Reg Strength 20","l");
    leg->AddEntry(NoneHist,"Reg Strength 0","l");
    leg->AddEntry(HistoData,"Data","l");
    leg->Draw();

    c->SaveAs("RegStr.png");
}

void PMNSParVal(){
    TFile* DataOsc = TFile::Open("EventRates/AllSampleTwentyTrueIndChanOsc.root");
    TFile* DataUnosc = TFile::Open("EventRates/AllSampleTwentyTrueIndChanUnosc.root");
    TFile* ColMap = TFile::Open("AllChanAllSamp(numu_to_numu)OscProbCol.root");

    std::vector<std::string> Samples = {"FHC_numu"};

    TFile* outFile = new TFile("RatioParam.root", "RECREATE");

    for (const auto& Sample : Samples) {     
        auto DataO = DataOsc->Get<TH1D>(("FD_" + Sample + "_FHC_numu_x_numu").c_str());
        auto DataU = DataUnosc->Get<TH1D>(("FD_" + Sample + "_FHC_numu_x_numu").c_str());

        TH1D* DataR = (TH1D*)DataO->Clone("DataRat");
        for (int i = 1; i <= DataO->GetNbinsX(); i++){
            double Oval = DataO->GetBinContent(i);
            double Uval = DataU->GetBinContent(i);
            if (Uval == 0) DataR->SetBinContent(i, 0);
            DataR->SetBinContent(i, Oval / Uval);
        }

        TCanvas *c = new TCanvas(("PostPred" + Sample).c_str(),"Overlay");
        ColMap->GetObject("(numu_to_numu)", c);
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