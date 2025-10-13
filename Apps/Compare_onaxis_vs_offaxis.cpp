#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"          // <-- needed for TMatrixD
#include <iostream>
#include <chrono>
#include "TLatex.h"
#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include "TPaveText.h"
#include <TColor.h>
#include <TMath.h>
#include <cmath> 
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <yaml-cpp/yaml.h>
#include <vector>
#include <ctime>
#include <TLine.h>
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"
#include <vector>
#include <string>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <random>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TKey.h"
#include "TCollection.h"

#include <iostream>
#include <string>
#include <regex>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>                 // <-- needed for std::set

#include <yaml-cpp/yaml.h>     // <-- needed for YAML::Node




int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <on plus off axis.root> <on axis.root>\n";
        return 1;
    }

    std::string fileA = argv[1];
    std::string fileB = argv[2];

    TFile* fA = TFile::Open(fileA.c_str(), "READ");
    TFile* fB = TFile::Open(fileB.c_str(), "READ");
    if (!fA || fA->IsZombie() || !fB || fB->IsZombie()) {
        std::cerr << "Error opening input files\n";
        return 1;
    }

    std::string pdfOut = "slice_ratios.pdf";
    TCanvas* c = new TCanvas("c","c",1200,2400);
    c->Print((pdfOut + "[").c_str());

    std::regex eventPattern("^h_event_slice_[0-9]+$");

    // loop over event histograms in file A
    TIter nextkey(fA->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        std::string name = key->GetName();
        name = name.substr(0, name.find(";")); // drop cycle suffix
        if (!std::regex_match(name, eventPattern)) continue;

        // matching posterior histogram
        std::string postName = std::regex_replace(name, std::regex("h_event"), "h_post");

        // get histograms from both files
        TH1D* h_event_A = (TH1D*) fA->Get(name.c_str());
        TH1D* h_event_B = (TH1D*) fB->Get(name.c_str());
        TH1D* h_post_A  = (TH1D*) fA->Get(postName.c_str());
        TH1D* h_post_B  = (TH1D*) fB->Get(postName.c_str());

        if (!h_event_A || !h_event_B || !h_post_A || !h_post_B) {
            std::cerr << "Skipping " << name << " (missing in one file)\n";
            continue;
        }

        // Get slice label from histogram title
        std::string sliceLabel = h_event_A->GetTitle();
        std::cout << "Hist name: " << h_event_A->GetName()
          << " | title: " << h_event_A->GetTitle()
          << " | x-axis title: " << h_event_A->GetXaxis()->GetTitle()
          << std::endl;

        // Determine x-axis range based on first and last bin of this slice
        //double xmin = h_event_A->GetBinLowEdge(1);
        //double xmax = h_event_A->GetBinLowEdge(h_event_A->GetNbinsX()) + h_event_A->GetBinWidth(h_event_A->GetNbinsX());

        // Determine x-axis range based on first and last *filled* bin in either file
        int firstFilled = -1;
        int lastFilled = -1;
        int nbins = h_event_A->GetNbinsX();

        

        // Look at event histograms A and B
        for (int i = 1; i <= nbins; ++i) {
            double vA = h_event_A->GetBinContent(i);
            double vB = (i <= h_event_B->GetNbinsX()) ? h_event_B->GetBinContent(i) : 0.0;
            if ((std::isfinite(vA) && vA != 0) || (std::isfinite(vB) && vB != 0)) {
                if (firstFilled < 0) firstFilled = i;
                lastFilled = i;
            }
        }

        // Also consider posterior histograms (in case event histograms are empty but posteriors are not)
        for (int i = 1; i <= h_post_A->GetNbinsX(); ++i) {
            double vA = h_post_A->GetBinContent(i);
            double vB = (i <= h_post_B->GetNbinsX()) ? h_post_B->GetBinContent(i) : 0.0;
            if ((std::isfinite(vA) && vA != 0) || (std::isfinite(vB) && vB != 0)) {
                if (firstFilled < 0) firstFilled = i;
                lastFilled = i;
            }
        }

        // Fallback: if no bins filled, use full range
        if (firstFilled < 0) {
            firstFilled = 1;
            lastFilled  = nbins;
        }

        // Now compute xmin/xmax from bin edges
        double xmin = h_event_A->GetXaxis()->GetBinLowEdge(firstFilled);
        double xmax = h_event_A->GetXaxis()->GetBinUpEdge(lastFilled);

        // --- compute y-axis range dynamically ---
        double ymin_ev = std::numeric_limits<double>::max();
        double ymax_ev = -std::numeric_limits<double>::max();

        for (int i = firstFilled; i <= lastFilled; ++i) {
            double vA = h_event_A->GetBinContent(i);
            double vB = (i <= h_event_B->GetNbinsX()) ? h_event_B->GetBinContent(i) : 0.0;
            if (std::isfinite(vA)) {
                ymin_ev = std::min(ymin_ev, vA);
                ymax_ev = std::max(ymax_ev, vA);
            }
            if (std::isfinite(vB)) {
                ymin_ev = std::min(ymin_ev, vB);
                ymax_ev = std::max(ymax_ev, vB);
            }
        }

        // fallback if histograms are empty
        if (ymax_ev <= -std::numeric_limits<double>::max()) {
            ymin_ev = 0.0;
            ymax_ev = 1.0;
        }

        // Ensure ymin is 0 or slightly below the smallest bin (for nice plots)
        if (ymin_ev > 0.0) ymin_ev = 0.0;

        // Add 10% headroom to ymax
        ymax_ev *= 1.1;

        h_event_A->GetYaxis()->SetRangeUser(ymin_ev, ymax_ev);
        h_event_B->GetYaxis()->SetRangeUser(ymin_ev, ymax_ev);

        // Clear canvas and divide into 3 pads
        c->Clear();
        c->Divide(1,4);

        // Add slice label as a title box at top of the canvas
        c->cd();
        TPaveText *titleBox = new TPaveText(0.1, 0.95, 0.9, 0.99, "NDC");
        titleBox->SetFillColor(0);
        titleBox->SetFillStyle(0);
        titleBox->SetBorderSize(0);
        titleBox->SetTextAlign(22); // centered
        titleBox->SetTextSize(0.04);
        titleBox->AddText(Form("Slice: %s", sliceLabel.c_str()));
        titleBox->Draw();


        // -------------------
        // Left pad: Event histograms
        // -------------------
        c->cd(1);
        h_event_A->SetLineColor(kRed);
        h_event_A->SetMarkerColor(kRed);
        h_event_A->SetMarkerStyle(20);

        h_event_B->SetLineColor(kBlue);
        h_event_B->SetMarkerColor(kBlue);
        h_event_B->SetMarkerStyle(21);

        h_event_A->GetXaxis()->SetRangeUser(xmin, xmax);
        h_event_B->GetXaxis()->SetRangeUser(xmin, xmax);

        h_post_A->GetXaxis()->SetRangeUser(xmin, xmax);
        h_post_B->GetXaxis()->SetRangeUser(xmin, xmax);

        h_event_A->SetTitle("");
        h_event_A->GetYaxis()->SetTitle("Counts");



        // Draw with error bars only (no connecting line)
        h_event_A->Draw("EZ");        
        h_event_B->Draw("EZ SAME");

        TLegend* leg1 = new TLegend(0.55,0.75,0.88,0.88);
        leg1->AddEntry(h_event_A, "On + Off axis", "lep");
        leg1->AddEntry(h_event_B, "On axis", "lep");
        leg1->Draw();
        c->cd(2);
        // Clone one histogram as a container for the ratio
        // Clone histogram to hold ratio
        TH1D* h_err_ratio = (TH1D*) h_event_A->Clone("h_err_ratio");
        h_err_ratio->Reset(); // clear contents but keep binning
        h_err_ratio->SetTitle("Ratio of Errors;X axis;#sigma_{A}/#sigma_{B}");

        for (int i = 1; i <= h_event_A->GetNbinsX(); i++) {
            double errA = h_event_A->GetBinError(i);
            double errB = h_event_B->GetBinError(i);

            double ratio = (errB > 0) ? errA / errB : 0.0;
            h_err_ratio->SetBinContent(i, ratio);
            h_err_ratio->SetBinError(i, 0.0); // no uncertainty on ratio
        }
        h_err_ratio->SetMaximum(0.55); 
        h_err_ratio->SetMinimum(0.65); 
        h_err_ratio->GetXaxis()->SetRangeUser(xmin, xmax);

        h_err_ratio->SetLineColor(kBlack);
        h_err_ratio->SetMarkerStyle(20);
        h_err_ratio->Draw("P");  // draw points only

                // -------------------
        // Middle pad: Posterior uncertainties
        // -------------------
        c->cd(3);
        TH1D* h_post_plot_A = (TH1D*) h_post_A->Clone("h_post_plot_A");
        TH1D* h_post_plot_B = (TH1D*) h_post_B->Clone("h_post_plot_B");

        h_post_plot_A->SetLineColor(kRed);   h_post_plot_A->SetMarkerColor(kRed); h_post_plot_A->SetMarkerStyle(20);
        h_post_plot_B->SetLineColor(kBlue);  h_post_plot_B->SetMarkerColor(kBlue); h_post_plot_B->SetMarkerStyle(21);

        h_post_plot_A->GetXaxis()->SetRangeUser(xmin, xmax);
        h_post_plot_B->GetXaxis()->SetRangeUser(xmin, xmax);

        h_post_plot_A->SetTitle("Posterior Uncertainties");
        h_post_plot_A->GetYaxis()->SetTitle("σ");
        h_post_plot_A->GetXaxis()->SetTitle(h_event_A->GetXaxis()->GetTitle());
        h_post_plot_A->GetXaxis()->SetRangeUser(xmin, xmax);
        h_post_plot_B->GetXaxis()->SetRangeUser(xmin, xmax);

        h_post_plot_A->Draw("E");
        h_post_plot_B->Draw("E SAME");

        TLegend* leg2 = new TLegend(0.55,0.75,0.88,0.88);
        leg2->AddEntry(h_post_plot_A, "On + Off axis", "lep");
        leg2->AddEntry(h_post_plot_B, "On axis", "lep");
        leg2->Draw();

        // -------------------
        // Right pad: Posterior uncertainty ratio (scatter)
        // -------------------
        /*
        c->cd(3);
        TH1D* h_unc_ratio = (TH1D*) h_post_A->Clone((postName + "_unc_ratio").c_str());
        h_unc_ratio->Reset();
        for (int b = 1; b <= h_post_A->GetNbinsX(); ++b) {
            double valA = h_post_A->GetBinError(b);
            double valB = h_post_B->GetBinError(b);
            double ratio = (valB > 0) ? valA / valB : 0.0;
            h_unc_ratio->SetBinContent(b, ratio);
        }

        h_unc_ratio->SetMarkerStyle(20);
        h_unc_ratio->SetMarkerColor(kMagenta-4);
        h_unc_ratio->SetTitle("Posterior Uncertainty Ratio");
        h_unc_ratio->GetYaxis()->SetTitle("On + Off Axis / On Axis");
        h_unc_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
        h_unc_ratio->Draw("P"); // scatter plot
        */
        // --- right pad: posterior uncertainty ratio as scatter ---
        // --- right pad: posterior uncertainty ratio as scatter ---
        c->cd(3);

        // Fill x_vals and y_vals
        std::vector<double> x_vals, y_vals;
        for (int b = 1; b <= h_post_A->GetNbinsX(); ++b) {
            double binCenter = h_post_A->GetBinCenter(b);
            double errA = h_post_A->GetBinError(b);
            double errB = h_post_B->GetBinError(b);
            double ratio = (errB > 0) ? errA / errB : 0.0;
            x_vals.push_back(binCenter);
            y_vals.push_back(ratio);
        }

        
        // Create scatter graph
        TGraph* g_ratio = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());
        g_ratio->SetMarkerStyle(21);
        g_ratio->SetMarkerColor(kBlack);
        g_ratio->SetLineColor(kBlack);
        g_ratio->SetTitle("Posterior Uncertainty Ratio");
        g_ratio->GetHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);
        gPad->Modified(); gPad->Update();
        g_ratio->GetXaxis()->SetTitle(h_post_A->GetXaxis()->GetTitle());
        g_ratio->GetYaxis()->SetTitle("On + Off Axis / On Axis");

        // Draw axes first
        g_ratio->Draw("AP");

        // Now get axis limits
        double xmin_canvas = g_ratio->GetXaxis()->GetXmin();
        double xmax_canvas = g_ratio->GetXaxis()->GetXmax();
        double ymin = g_ratio->GetYaxis()->GetXmin();
        double ymax = g_ratio->GetYaxis()->GetXmax();

        // Draw pale green/red boxes in **graph coordinates**
        TBox* greenBox = new TBox(xmin, ymin, xmax, 1.0);
        greenBox->SetFillColorAlpha(kGreen+2, 0.2);
        greenBox->SetLineStyle(0);
        greenBox->Draw("SAME");

        TBox* redBox = new TBox(xmin, 1.0, xmax, ymax);
        redBox->SetFillColorAlpha(kRed+2, 0.2);
        redBox->SetLineStyle(0);
        redBox->Draw("SAME");

        // Draw horizontal dotted line at y=1
        TLine* line = new TLine(xmin, 1.0, xmax, 1.0);
        line->SetLineColor(kBlack);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("SAME");

                // --- pad 4: posterior uncertainty ratio as scatter ---
        c->cd(4);

        

        TGraph*  g_ratio2 = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());
        g_ratio2->SetTitle("Posterior Uncertainty Ratio;X;#sigma_{A}/#sigma_{B}");
        g_ratio2->SetMarkerStyle(20);
        g_ratio2->SetMarkerColor(kMagenta+2);
        g_ratio2->GetHistogram()->GetXaxis()->SetRangeUser(xmin, xmax);
        gPad->Modified(); gPad->Update();
        g_ratio2->Draw("AP");

        // Add slice label at top
        c->cd();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextAlign(22);
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.5, 0.96, sliceLabel.c_str());

        // Save this page
        c->Print(pdfOut.c_str());

        // cleanup
        delete leg1;
        delete leg2;
        delete h_post_plot_A;
        delete h_post_plot_B;
        delete g_ratio;
    }

    c->Print((pdfOut + "]").c_str());

    fA->Close();
    fB->Close();
    delete c;
    return 0;
}
 