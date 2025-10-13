#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " file1.root file2.root" << std::endl;
        return 1;
    }

    std::string fileA = argv[1];
    std::string fileB = argv[2];

    TFile *fA = TFile::Open(fileA.c_str());
    TFile *fB = TFile::Open(fileB.c_str());
    if (!fA || !fB || fA->IsZombie() || fB->IsZombie()) {
        std::cerr << "Error: Could not open input files!" << std::endl;
        return 1;
    }

    // Collect event slice names automatically from fileA
    std::vector<std::string> sliceNames;
    TIter nextKey(fA->GetListOfKeys());
    while (TKey *key = (TKey*) nextKey()) {
        if (std::string(key->GetClassName()) == "TH1D") {
            std::string name = key->GetName();
            if (name.find("h_event_slice_") == 0) {
                if (fB->Get(name.c_str())) { // only keep if present in fileB
                    sliceNames.push_back(name);
                }
            }
        }
    }

    std::cout << "Found " << sliceNames.size() << " common slices to compare." << std::endl;
    std::string outFile = "all_slices_comparison.pdf";
    bool firstPage = true;

    for (const auto &evName : sliceNames) {
        std::string idx = evName.substr(std::string("h_event_slice_").size());
        std::string postName = "h_post_slice_" + idx;

        
        


        TH1D *h_event_A = (TH1D*) fA->Get(evName.c_str());
        TH1D *h_event_B = (TH1D*) fB->Get(evName.c_str());
        TH1D *h_post_A  = (TH1D*) fA->Get(postName.c_str());
        TH1D *h_post_B  = (TH1D*) fB->Get(postName.c_str());

        if (!h_event_A || !h_event_B || !h_post_A || !h_post_B) {
            std::cerr << "Skipping " << evName << " (missing post hist in one file)" << std::endl;
            continue;
        }

        // ✅ Normalize histograms
        double intA = h_event_A->Integral();
        double intB = h_event_B->Integral();
        if (intA > 0) h_event_A->Scale(1.0 / intA);
        if (intB > 0) h_event_B->Scale(1.0 / intB);

        double postIntA = h_post_A->Integral();
        double postIntB = h_post_B->Integral();
        if (postIntA > 0) h_post_A->Scale(1.0 / postIntA);
        if (postIntB > 0) h_post_B->Scale(1.0 / postIntB);

        // ✅ Find first and last filled bins
        int nbins = h_post_A->GetNbinsX();
        int firstFilled = 0, lastFilled = 0;
        double threshold = 1e-12;

        for (int b = 1; b <= nbins; ++b) {
            double sum = h_post_A->GetBinContent(b) + h_post_B->GetBinContent(b)
                       + h_event_A->GetBinContent(b) + h_event_B->GetBinContent(b);
            if (sum > threshold) { firstFilled = b; break; }
        }
        for (int b = nbins; b >= 1; --b) {
            double sum = h_post_A->GetBinContent(b) + h_post_B->GetBinContent(b)
                       + h_event_A->GetBinContent(b) + h_event_B->GetBinContent(b);
            if (sum > threshold) { lastFilled = b; break; }
        }

        double xlow  = h_post_A->GetXaxis()->GetBinLowEdge(firstFilled > 0 ? firstFilled : 1);
        double xhigh = h_post_A->GetXaxis()->GetBinUpEdge(lastFilled > 0 ? lastFilled : nbins);
        h_event_A->GetXaxis()->SetRangeUser(xlow, xhigh);
        h_event_B->GetXaxis()->SetRangeUser(xlow, xhigh);
        h_post_A->GetXaxis()->SetRangeUser(xlow, xhigh);
        h_post_B->GetXaxis()->SetRangeUser(xlow, xhigh);

        // ✅ Make ratio histogram
        TH1D *h_err_ratio = (TH1D*) h_event_A->Clone((evName + "_errRatio").c_str());
        h_err_ratio->Reset("ICES");
        double minRatio = 1e9, maxRatio = -1e9;
        for (int b = 1; b <= nbins; ++b) {
            double errA = h_event_A->GetBinError(b);
            double errB = h_event_B->GetBinError(b);
            double ratio = (errB > 0) ? errA / errB : 0.0;
            h_err_ratio->SetBinContent(b, ratio);
            if (ratio > 0) {
                minRatio = std::min(minRatio, ratio);
                maxRatio = std::max(maxRatio, ratio);
            }
        }
        if (minRatio < maxRatio) {
            h_err_ratio->SetMinimum(0.9 * minRatio);
            h_err_ratio->SetMaximum(1.1 * maxRatio);
        }

        // ✅ Draw
        TCanvas *c = new TCanvas((evName + "_canvas").c_str(), evName.c_str(), 1000, 800);
        c->Divide(1, 3);

        c->cd(1);
        h_event_A->SetLineColor(kBlue);
        h_event_B->SetLineColor(kRed);
        h_event_A->Draw("HIST E");
        h_event_B->Draw("HIST E SAME");

        auto leg1 = new TLegend(0.7, 0.75, 0.9, 0.9);
        leg1->AddEntry(h_event_A, "File 1", "l");
        leg1->AddEntry(h_event_B, "File 2", "l");
        leg1->Draw();

        c->cd(2);
        h_post_A->SetLineColor(kBlue);
        h_post_B->SetLineColor(kRed);
        h_post_A->Draw("HIST E");
        h_post_B->Draw("HIST E SAME");

        auto leg2 = new TLegend(0.7, 0.75, 0.9, 0.9);
        leg2->AddEntry(h_post_A, "Post File 1", "l");
        leg2->AddEntry(h_post_B, "Post File 2", "l");
        leg2->Draw();

        c->cd(3);
        h_err_ratio->SetLineColor(kBlack);
        h_err_ratio->SetTitle("Error Ratio (File1/File2)");
        h_err_ratio->Draw("HIST");

        if (firstPage) {
            c->SaveAs((outFile + "(").c_str()); // open multi-page PDF
            firstPage = false;
        } else {
            c->SaveAs(outFile.c_str()); // add a new page
        }

        //c->SaveAs((evName + "_comparison.pdf").c_str());

        delete h_err_ratio;
    }

    if (!firstPage) {
    TCanvas dummy;
    dummy.SaveAs((outFile + ")").c_str()); // close the PDF
    }


    fA->Close();
    fB->Close();
    return 0;
}
