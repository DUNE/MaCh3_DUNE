#include <TFile.h>
#include <TSpline.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TKey.h>
#include <TLatex.h>
#include <TCollection.h>
#include <TAxis.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

void plotSplines(const std::string& inputFile = "OffAxisND_splines.root", const std::string& outputPDF = "splines.pdf") {

    TFile* f = TFile::Open(inputFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open file: " << inputFile << std::endl;
        return;
    }

    std::vector<std::string> splineNames = {
        "maqe", "vecffqeshape", "pauli", "mares", "mvres", "mancres", "mvncres",
        "thetadelta", "ahtby", "bhtby", "cv1uby", "cv2uby", "frcexpi", "frelaspi",
        "frinelpi", "frabspi", "frpiprodpi", "frcexn", "frelasn", "frineln",
        "frabsn", "frpiprodn", "e2anu", "e2bnu", "e2anubar", "e2bnubar",
        "c12nu", "c12nubar", "nuncc2", "nuncc3", "nupcc2", "nupcc3", "nunpcc1",
        "nunnc1", "nunnc2", "nunnc3", "nupnc1", "nupnc2", "nupnc3",
        "nubarncc1", "nubarncc2", "nubarncc3", "nubarpcc1", "nubarpcc2", "nubarpcc3",
        "nubarnnc1", "nubarnnc2", "nubarnnc3", "nubarpnc1", "nubarpnc2", "nubarpnc3",
        "berpaa", "berpab", "berpad", "nuexsec", "nuemuxsec"
    };

    std::vector<int> colors = {
        kBlue, kRed, kGreen+2, kMagenta+1, kOrange+7, kCyan+1,
        kViolet+5, kSpring-1, kAzure+1, kPink+6, kTeal+3, kYellow+2
    };

    const int    npoints   = 100;
    const double tolerance = 1e-6;
    const double SigmaKnots[] = {-3, -2, -1, 0, 1, 2, 3};

    TCanvas* c = new TCanvas("c", "c", 1200, 600);
    c->Print((outputPDF + "[").c_str());

    for (const auto& sysname : splineNames) {

        // Collect all splines whose names contain this spline name
        std::vector<TSpline3*> splines;
        TIter next(f->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            std::string objName = key->GetName();
            if (objName.find(sysname) != std::string::npos) {
                TSpline3* sp = dynamic_cast<TSpline3*>(f->Get(objName.c_str()));
                if (sp) splines.push_back(sp);
            }
        }

        if (splines.empty()) continue;

        c->Clear();
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.05);
        c->SetTopMargin(0.12); // room for title

        // Build the title from the names of all splines on this page
        std::string titleStr = sysname; // + ":  ";
        // for (size_t i = 0; i < splines.size(); i++) {
        //     titleStr += splines[i]->GetName();
        //     if (i + 1 < splines.size()) titleStr += ",  ";
        // }

        std::vector<TGraph*> drawnGraphs;

        // First pass: compute global Y min/max across all drawn graphs
        double ymin =  1e99;
        double ymax = -1e99;

        for (auto* sp : splines) {
            double nominal = sp->Eval(0.0);
            if (nominal == 0.) continue;
            bool all_one = true;
            for (int i = 0; i < npoints; i++) {
                double sigma = -3.0 + 6.0*i/(npoints-1);
                double ratio = sp->Eval(sigma) / nominal;
                if (std::fabs(ratio - 1.0) > tolerance) all_one = false;
                ymin = std::min(ymin, ratio);
                ymax = std::max(ymax, ratio);
            }
            if (all_one) continue; // still update ymin/ymax above but skip drawing below
        }

        // Add 10% padding
        double ypad = 0.10 * (ymax - ymin);
        ymin -= ypad;
        ymax += ypad;

        // Second pass: draw
        for (auto* sp : splines) {
            double nominal = sp->Eval(0.0);
            if (nominal == 0.) continue;

            bool all_one = true;
            for (int i = 0; i < npoints; i++) {
                double sigma = -3.0 + 6.0*i/(npoints-1);
                double ratio = sp->Eval(sigma) / nominal;
                if (std::fabs(ratio - 1.0) > tolerance) {
                    all_one = false;
                    break;
                }
            }
            if (all_one) continue;

            TGraph* g = new TGraph();
            for (int i = 0; i < npoints; i++) {
                double sigma = -3.0 + 6.0*i/(npoints-1);
                double ratio = sp->Eval(sigma) / nominal;
                g->SetPoint(i, sigma, ratio);
            }

            int color = colors[drawnGraphs.size() % colors.size()];
            g->SetLineWidth(2);
            g->SetLineColor(color);

            if (drawnGraphs.empty()) {
                g->SetTitle(";#sigma;Variation / Nominal");
                g->GetXaxis()->SetRangeUser(-3.5, 3.5);
                g->GetYaxis()->SetRangeUser(ymin, ymax);
                g->Draw("AL");
            } else {
                g->Draw("L SAME");
            }

            // Draw knot markers
            for (int k = 0; k < 7; k++) {
                double sigma = SigmaKnots[k];
                double ratio = sp->Eval(sigma) / nominal;
                TMarker* m = new TMarker(sigma, ratio, 20);
                m->SetMarkerColor(color);
                m->SetMarkerSize(0.8);
                m->Draw("same");
            }

            drawnGraphs.push_back(g);
        }

        if (!drawnGraphs.empty()) {
            // Draw title using TLatex in NDC
            TLatex title;
            title.SetNDC();
            title.SetTextSize(0.035);
            title.SetTextFont(62);
            title.SetTextAlign(22);
            title.DrawLatex(0.5, 0.97, titleStr.c_str());

            c->Update();
            c->Print(outputPDF.c_str());
        }

        for (auto* g : drawnGraphs) delete g;
        drawnGraphs.clear();
    }

    c->Print((outputPDF + "]").c_str());

    f->Close();
    delete f;
    delete c;

    std::cout << "Saved: " << outputPDF << std::endl;
}