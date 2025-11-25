// #include <TFile.h>
// #include <TH1D.h>
// #include <TCanvas.h>
// #include <TLegend.h>
// #include <TGraphErrors.h>
// #include <TKey.h>
// #include <TDirectory.h>
// #include <TString.h>
// #include <iostream>
// #include <vector>
// #include <algorithm>
// #include <cmath>
// #include <functional>


// auto addStatErrors = [](TH1D* meanHist, TH1D* sigmaHist, const std::vector<TH1D*>& samples){
//     if (!meanHist || !sigmaHist || samples.empty()) return (TH1D*)nullptr;

//     // Clone the posterior sigma before modifying it, so we can draw both later
//     TH1D* hPosteriorOnly = (TH1D*)sigmaHist->Clone(TString(sigmaHist->GetName()) + "_posteriorOnly");

//     int nb = meanHist->GetNbinsX();
//     for (int b=1; b<=nb; ++b){
//         double postSigma = sigmaHist->GetBinContent(b);

//         double sumStat2 = 0;
//         for (auto* h : samples) {
//             double statErr = h->GetBinError(b);
//             sumStat2 += statErr * statErr;
//         }
//         double meanStatErr = sqrt(sumStat2 / samples.size());
//         double totalSigma = sqrt(postSigma*postSigma + meanStatErr*meanStatErr);

//         sigmaHist->SetBinContent(b, totalSigma);
//     }

//     return hPosteriorOnly; // return a clone of the original sigma
// };



// void fd_nd_postpred(const char* filename="/scratch/abipeake/October_chains/EventRates_1D_adapted_70templateparams_25erecobins/posteriorpredictive.root"){ //scratch/abipeake/October_chains/EventRates_1D_adapted_70templateparams_25erecobins/posteriorpredictive.root") {
//     TFile* f = TFile::Open(filename);
//     if (!f || f->IsZombie()) {
//         std::cerr << "❌ Could not open file " << filename << "\n";
//         return;
//     }

//     // Helper to find histograms in a directory with prefix + tag
//     auto getHists = [&](const TString& tag, const TString& prefix, const TString& dirHint) {
//         std::vector<TH1D*> vec;

//         std::function<void(TDirectory*)> recurse = [&](TDirectory* dir) {
//             TIter next(dir->GetListOfKeys());
//             TKey* key;
//             while ((key = (TKey*)next())) {
//                 TObject* obj = key->ReadObj();
//                 if (obj->InheritsFrom(TDirectory::Class())) {
//                     recurse((TDirectory*)obj);
//                 } else if (obj->InheritsFrom(TH1D::Class())) {
//                     TString name = obj->GetName();
//                     TString path = dir->GetPath();
//                     if (name.BeginsWith(prefix) && name.Contains(tag) && path.Contains(dirHint)) {
//                         vec.push_back((TH1D*)obj);
//                     }
//                 }
//             }
//         };

//         recurse(f);

//         std::sort(vec.begin(), vec.end(),
//                   [](TH1D* a, TH1D* b){ return TString(a->GetName()) < TString(b->GetName()); });

//         if (!vec.empty()) {
//             TString dirName = vec.front()->GetDirectory()->GetName();
//             std::cout << "✅ Found " << vec.size()
//                       << " histograms with prefix \"" << prefix
//                       << "\" and tag \"" << tag
//                       << "\" in directory \"" << dirName << "\"\n";
//         } else {
//             std::cout << "⚠️  No histograms found with prefix \"" << prefix
//                       << "\" and tag \"" << tag
//                       << "\" in directory \"" << dirHint << "\"\n";
//         }

//         return vec;
//     };

//     // --- compute mean and sigma histograms ---
//     auto makeMeanSigma = [](const std::vector<TH1D*>& hists, const char* base){
//         if (hists.empty()) return std::pair<TH1D*,TH1D*>(nullptr, nullptr);
//         TH1D* hMean = (TH1D*)hists[0]->Clone(TString(base)+"_mean");
//         TH1D* hSigma = (TH1D*)hists[0]->Clone(TString(base)+"_sigma");
//         hMean->Reset(); hSigma->Reset();
//         int nb = hMean->GetNbinsX();
//         for (int b=1; b<=nb; ++b) {
//             double sum=0, sum2=0;
//             for (auto* h : hists) {
//                 double v = h->GetBinContent(b);
//                 sum += v; sum2 += v*v;
//             }
//             double mean = sum / hists.size();
//             double var = sum2 / hists.size() - mean*mean;
//             hMean->SetBinContent(b, mean);
//             hSigma->SetBinContent(b, var>0 ? sqrt(var) : 0);
//         }
//         return std::make_pair(hMean, hSigma);
//     };

//     // --- divide by bin width ---
//     auto divideByBinWidth = [](TH1D* h){
//         if (!h) return;
//         for (int b=1; b<=h->GetNbinsX(); ++b) {
//             double width = h->GetBinWidth(b);
//             if (width > 0) h->SetBinContent(b, h->GetBinContent(b)/width);
//         }
//     };

//     // --- find histograms ---
//     auto NDhists = getHists("", "ND_FHC_CCnumu_posterior_sample", "ND");
//     auto FDhists = getHists("", "FHC_numu_posterior_sample_", "Other"); // FD histos are in Unknown directory

    
//     auto [NDmean, NDsigma] = makeMeanSigma(NDhists, "ND");
//     auto [FDmean, FDsigma] = makeMeanSigma(FDhists, "FD");

//     //TH1D* NDposteriorOnly = addStatErrors(NDmean, NDsigma, NDhists);
//     //TH1D* FDposteriorOnly = addStatErrors(FDmean, FDsigma, FDhists);
//     TH1D* NDposteriorOnly = NDsigma;  // reuse sigma as-is
//     TH1D* FDposteriorOnly = FDsigma;


//     divideByBinWidth(NDmean);
//     divideByBinWidth(NDsigma);
//     divideByBinWidth(FDmean);
//     divideByBinWidth(FDsigma);

//     if (!NDmean && !FDmean) {
//         std::cerr << "❌ No ND or FD histograms found!\n";
//         return;
//     }

//     // --- ND canvas ---
//     if (NDmean) {
//         TCanvas* cND = new TCanvas("cND","Posterior Predictive ND",800,600);
//         cND->SetGrid();

//         NDmean->SetLineColor(kBlue+1);
//         NDmean->SetLineWidth(2);
//         NDmean->SetTitle("Posterior Predictive ND;Reco Neutrino Energy [GeV];Events / GeV");
//         NDmean->Draw("hist");

//         if (NDsigma) {
//             // Light blue: posterior-only band
//             if (NDposteriorOnly) {
//                 TGraphErrors* grNDpost = new TGraphErrors(NDmean->GetNbinsX());
//                 for (int b=1; b<=NDmean->GetNbinsX(); ++b){
//                     double x = NDmean->GetBinCenter(b);
//                     double y = NDmean->GetBinContent(b);
//                     double e = NDposteriorOnly->GetBinContent(b);
//                     grNDpost->SetPoint(b-1, x, y);
//                     grNDpost->SetPointError(b-1, 0, e);
//                 }
//                 grNDpost->SetFillColorAlpha(kBlue-9, 0.2);
//                 grNDpost->Draw("E3 same");
//             }

//             // Darker blue: total uncertainty band
//             TGraphErrors* grNDtotal = new TGraphErrors(NDmean->GetNbinsX());
//             for (int b=1; b<=NDmean->GetNbinsX(); ++b){
//                 double x = NDmean->GetBinCenter(b);
//                 double y = NDmean->GetBinContent(b);
//                 double e = NDsigma->GetBinContent(b);
//                 grNDtotal->SetPoint(b-1, x, y);
//                 grNDtotal->SetPointError(b-1, 0, e);
//             }
//             grNDtotal->SetFillColorAlpha(kBlue-3, 0.35);
//             grNDtotal->Draw("E3 same");
//         }


//         TLegend* legND = new TLegend(0.62,0.73,0.88,0.88);
//         legND->AddEntry(NDmean,"ND mean","l");
//         //legND->AddEntry((TObject*)0,"Shaded bands:","");
//         legND->AddEntry((TObject*)0," posterior","");
//         //legND->AddEntry((TObject*)0,"dark = total (posterior + stat)","");
//         legND->Draw();


//         cND->SaveAs("ND_posterior_70templateparams_25enurecndobservables.pdf");
//         std::cout << "→ Saved ND plot: ND_posterior_70templateparams_25enurecndobservables.pdf\n";
//     }

//     // --- FD canvas ---
//     if (FDmean) {
//         TCanvas* cFD = new TCanvas("cFD","Posterior Predictive FD",800,600);
//         cFD->SetGrid();

//         FDmean->SetLineColor(kRed+1);
//         FDmean->SetLineWidth(2);
//         FDmean->SetTitle("Posterior Predictive FD;Reco Neutrino Energy [GeV];Events / GeV");
//         FDmean->Draw("hist");

//         if (FDsigma) {
//             TGraphErrors* grFD = new TGraphErrors(FDmean->GetNbinsX());
//             for (int b=1;b<=FDmean->GetNbinsX();++b){
//                 double x=FDmean->GetBinCenter(b);
//                 double y=FDmean->GetBinContent(b);
//                 double e=FDsigma->GetBinContent(b);
//                 grFD->SetPoint(b-1,x,y);
//                 grFD->SetPointError(b-1,0,e);
//             }
//             grFD->SetFillColorAlpha(kRed-9, 0.3);
//             grFD->Draw("E1 same");
//         }

//         TLegend* legFD = new TLegend(0.65,0.75,0.88,0.88);
//         legFD->AddEntry(FDmean,"FD mean #pm 1#sigma","lep");
//         legFD->Draw();


//         // --- Save output ROOT file ---
//         TFile* fout = new TFile("posterior_predictive_outputtuesday.root", "RECREATE");

//         // Save ND results
//         if (NDmean) NDmean->Write();
//         if (NDsigma) NDsigma->Write();
//         if (NDposteriorOnly) NDposteriorOnly->Write();

//         // Save FD results
//         if (FDmean) FDmean->Write();
//         if (FDsigma) FDsigma->Write();
//         if (FDposteriorOnly) FDposteriorOnly->Write();

//         // Optionally save the canvases and graphs too
//         if (gROOT->FindObject("cND")) ((TCanvas*)gROOT->FindObject("cND"))->Write();
//         if (gROOT->FindObject("cFD")) ((TCanvas*)gROOT->FindObject("cFD"))->Write();

//         fout->Close();
//         std::cout << "💾 Saved all histograms and canvases to posterior_predictive_output.root\n";


//         cFD->SaveAs("FD_posterior_70templateparams_25enurecndobservables.pdf");
//         std::cout << "→ Saved FD plot: FD_posterior_70templateparams_25enurecndobservables.pdf\n";
//     }

//     std::cout << "✅ Finished processing posterior predictive plots.\n";
// }


#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <TStyle.h>
#include <TPad.h>

void fd_nd_postpred(const char* filename="/scratch/abipeake/October_chains/EventRates_EnuEnubias_adapting_biggerstepscale/posteriorpred.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.045);

    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "❌ Could not open file " << filename << "\n";
        return;
    }

    // Helper to find histograms recursively
    auto getHists = [&](const TString& tag, const TString& prefix, const TString& dirHint) {
        std::vector<TH1D*> vec;

        std::function<void(TDirectory*)> recurse = [&](TDirectory* dir) {
            if (!dir) return;
            TIter next(dir->GetListOfKeys());
            TKey* key;
            while ((key = (TKey*)next())) {
                TObject* obj = key->ReadObj();
                if (obj->InheritsFrom(TDirectory::Class())) {
                    recurse((TDirectory*)obj);
                } else if (obj->InheritsFrom(TH1D::Class())) {
                    TString name = obj->GetName();
                    TString path = dir->GetPath();
                    if (name.BeginsWith(prefix) && name.Contains(tag) && path.Contains(dirHint)) {
                        TH1D* h = (TH1D*)obj;
                        h->SetDirectory(nullptr);
                        vec.push_back(h);
                    }
                }
            }
        };

        recurse(f);

        std::sort(vec.begin(), vec.end(),
                  [](TH1D* a, TH1D* b){ return TString(a->GetName()) < TString(b->GetName()); });

        if (!vec.empty()) {
            std::cout << "✅ Found " << vec.size()
                      << " histograms with prefix \"" << prefix
                      << "\" and tag \"" << tag
                      << "\" in directory hint \"" << dirHint << "\"\n";
        } else {
            std::cout << "⚠️  No histograms found with prefix \"" << prefix
                      << "\" and tag \"" << tag
                      << "\" in directory hint \"" << dirHint << "\"\n";
        }

        return vec;
    };

    // Mean + Sigma
    auto makeMeanSigma = [](const std::vector<TH1D*>& hists, const char* base){
        if (hists.empty()) return std::pair<TH1D*,TH1D*>(nullptr, nullptr);
        TH1D* hMean = (TH1D*)hists[0]->Clone(TString(base)+"_mean");
        TH1D* hSigma = (TH1D*)hists[0]->Clone(TString(base)+"_sigma");
        hMean->Reset(); hSigma->Reset();
        int nb = hMean->GetNbinsX();
        for (int b=1; b<=nb; ++b) {
            double sum=0, sum2=0;
            for (auto* h : hists) {
                double v = h->GetBinContent(b);
                sum += v; sum2 += v*v;
            }
            double mean = sum / hists.size();
            double var = sum2 / hists.size() - mean*mean;
            hMean->SetBinContent(b, mean);
            hSigma->SetBinContent(b, var>0 ? sqrt(var) : 0);
        }
        return std::make_pair(hMean, hSigma);
    };

    // Divide histograms by bin width
    auto divideByBinWidth = [](TH1D* h){
        if (!h) return;
        for (int b=1; b<=h->GetNbinsX(); ++b){
            double width = h->GetBinWidth(b);
            if (width>0) h->SetBinContent(b, h->GetBinContent(b)/width);
        }
    };

    auto makeFracPosteriorErr = [](TH1D* meanHist, TH1D* sigmaHist){
        if (!meanHist || !sigmaHist) return (TH1D*)nullptr;
        TH1D* hFrac = (TH1D*)meanHist->Clone(TString(meanHist->GetName()) + "_fracPosteriorErr");
        hFrac->Reset();
        hFrac->SetTitle("Fractional Posterior Uncertainty;Reco Neutrino Energy [GeV];Uncertainty [%]");

        int nb = meanHist->GetNbinsX();
        for (int b=1; b<=nb; ++b){
            double meanVal = meanHist->GetBinContent(b);
            double postErr = sigmaHist->GetBinContent(b); // <-- take the posterior error
            double fracErr = (meanVal != 0) ? 100.0 * postErr / meanVal : 0.0;
            hFrac->SetBinContent(b, fracErr);
        }
        return hFrac;
};


    // --- Get histograms ---
    auto NDhists = getHists("", "ND_FHC_CCnumu_posterior_sample", "ND");
    auto FDhists = getHists("", "FHC_numu_posterior_sample_", "Other");

    auto [NDmean, NDsigma] = makeMeanSigma(NDhists, "ND");
    auto [FDmean, FDsigma] = makeMeanSigma(FDhists, "FD");

    divideByBinWidth(NDmean);
    divideByBinWidth(NDsigma);
    divideByBinWidth(FDmean);
    divideByBinWidth(FDsigma);

    //TH1D* NDfracPosterior = makeFracPosteriorErr(NDmean, NDhists);
    //TH1D* FDfracPosterior = makeFracPosteriorErr(FDmean, FDhists);
    TH1D* NDfracPosterior = makeFracPosteriorErr(NDmean, NDsigma);
    TH1D* FDfracPosterior = makeFracPosteriorErr(FDmean, FDsigma);


    if (!NDmean && !FDmean) {
        std::cerr << "❌ No ND or FD histograms found!\n";
        f->Close();
        return;
    }

    // ---------------------- ND CANVAS ----------------------
    if (NDmean) {
        TCanvas* cND = new TCanvas("cND","Posterior Predictive ND",800,700);
        cND->Divide(1,2);

        // --- Top pad: posterior predictive ---
        TPad* pad1 = (TPad*)cND->cd(1);
        pad1->SetPad(0,0.3,1,1);
        pad1->SetBottomMargin(0.02);
        pad1->SetGrid();

        NDmean->SetLineColor(kBlue+1);
        NDmean->SetLineWidth(2);
        NDmean->SetTitle("Posterior Predictive ND;Reco Neutrino Energy [GeV];Events / GeV");
        NDmean->Draw("hist");

        if (NDsigma) {
            TGraphErrors* grND = new TGraphErrors(NDmean->GetNbinsX());
            for (int b=1;b<=NDmean->GetNbinsX();++b){
                double x=NDmean->GetBinCenter(b);
                double y=NDmean->GetBinContent(b);
                double e=NDsigma->GetBinContent(b);
                grND->SetPoint(b-1,x,y);
                grND->SetPointError(b-1,0,e);
            }
            grND->SetFillColorAlpha(kBlue-9, 0.3);
            grND->Draw("E3 same");
        }

        TLegend* legND = new TLegend(0.62,0.73,0.88,0.88);
        legND->AddEntry(NDmean,"ND mean #pm posterior #sigma","l");
        legND->Draw();

        // --- Bottom pad: fractional statistical error ---
        cND->cd(2);
        TPad* pad2 = (TPad*)gPad;
        pad2->SetPad(0,0,1,0.3);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.25);
        pad2->SetGrid();

        if (NDfracPosterior) {
            NDfracPosterior->SetLineColor(kBlue+1);
            NDfracPosterior->SetLineWidth(2);
            NDfracPosterior->SetMarkerStyle(20);
            NDfracPosterior->SetMarkerSize(0.7);
            NDfracPosterior->GetYaxis()->SetTitleSize(0.09);
            NDfracPosterior->GetYaxis()->SetLabelSize(0.08);
            NDfracPosterior->GetXaxis()->SetTitleSize(0.09);
            NDfracPosterior->GetXaxis()->SetLabelSize(0.08);
            NDfracPosterior->Draw("E1");
        }


        cND->SaveAs("ND_posterior_with_adaptedenubiaschain.pdf");
        std::cout << "→ Saved ND plot: ND_posterior_with_adaptedenubiaschain.pdf\n";
    }

    // ---------------------- FD CANVAS ----------------------
    if (FDmean) {
        TCanvas* cFD = new TCanvas("cFD","Posterior Predictive FD",800,700);
        cFD->Divide(1,2);

        // --- Top pad ---
        TPad* pad1 = (TPad*)cFD->cd(1);
        pad1->SetPad(0,0.3,1,1);
        pad1->SetBottomMargin(0.02);
        pad1->SetGrid();

        FDmean->SetLineColor(kRed+1);
        FDmean->SetLineWidth(2);
        FDmean->SetTitle("Posterior Predictive FD;Reco Neutrino Energy [GeV];Events / GeV");
        FDmean->Draw("hist");

        if (FDsigma) {
            TGraphErrors* grFD = new TGraphErrors(FDmean->GetNbinsX());
            for (int b=1;b<=FDmean->GetNbinsX();++b){
                double x=FDmean->GetBinCenter(b);
                double y=FDmean->GetBinContent(b);
                double e=FDsigma->GetBinContent(b);
                grFD->SetPoint(b-1,x,y);
                grFD->SetPointError(b-1,0,e);
            }
            grFD->SetFillColorAlpha(kRed-9, 0.3);
            grFD->Draw("E1 same");
        }

        TLegend* legFD = new TLegend(0.65,0.75,0.88,0.88);
        legFD->AddEntry(FDmean,"FD mean #pm posterior #sigma","l");
        legFD->Draw();

        // --- Bottom pad: fractional stat error ---
        cFD->cd(2);
        TPad* pad2 = (TPad*)gPad;
        pad2->SetPad(0,0,1,0.3);
        pad2->SetTopMargin(0.05);
        pad2->SetBottomMargin(0.25);
        pad2->SetGrid();

        if (FDfracPosterior) {
            FDfracPosterior->SetLineColor(kRed+1);
            FDfracPosterior->SetLineWidth(2);
            FDfracPosterior->SetMarkerStyle(20);
            FDfracPosterior->SetMarkerSize(0.7);
            FDfracPosterior->GetYaxis()->SetTitleSize(0.09);
            FDfracPosterior->GetYaxis()->SetLabelSize(0.08);
            FDfracPosterior->GetXaxis()->SetTitleSize(0.09);
            FDfracPosterior->GetXaxis()->SetLabelSize(0.08);
            FDfracPosterior->GetYaxis()->SetRangeUser(0,5);
            FDfracPosterior->Draw("");
        }


        cFD->SaveAs("FD_posterior_with_adaptedenubiaschain.pdf");
        std::cout << "→ Saved FD plot: FD_posterior_with_adaptedenubiaschain.pdf\n";
    }

    // --- Save output ROOT file ---
    TFile* fout = new TFile("posterior_predictive_output_withFracStat.root", "RECREATE");
    if (NDmean) NDmean->Write();
    if (NDsigma) NDsigma->Write();
    if (NDfracPosterior) NDfracPosterior->Write();
    if (FDmean) FDmean->Write();
    if (FDsigma) FDsigma->Write();
    if (FDfracPosterior) FDfracPosterior->Write();
    fout->Close();

    f->Close();

    std::cout << "💾 Saved all histograms and canvases to posterior_predictive_output_withFracStat.root\n";
    std::cout << "✅ Finished processing posterior predictive plots.\n";
}
