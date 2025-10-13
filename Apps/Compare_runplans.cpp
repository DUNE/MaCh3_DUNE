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
#include "TMatrixD.h"
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <yaml-cpp/yaml.h>
#include "manager/manager.h"

/*
std::vector<std::string> files_A = {
        "/scratch/abipeake/Off_axis_chains_new/RunPlan_10_90_afteradapted/RunPlan_10_90_afteradapted_chain_0_job_0_diagnostics.root",
        "/scratch/abipeake/Off_axis_chains_new/RunPlan_25_75_afteradapted_new/RunPlan_25_75_afteradapted_new_chain_0_job_0_diagnostics.root",
        "/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2_chain_0_job_0backup_diagnostics.root",
        "/scratch/abipeake/Off_axis_chains_new/RunPlan_75_25_afteradapted/RunPlan_75_25_afteradapted_chain_0_job_0_diagnostics.root",
        "/scratch/abipeake/Off_axis_chains_new/RunPlan_90_10_afteradapted/RunPlan_90_10_afteradapted_chain_0_job_0_diagnostics.root"
    };
    std::vector<std::string> files_B = {
        "/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/Onaxis_EventRates_Elep_erec_afteradaptive_chain_1_job_0test_diagnostics.root",
        "/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/Onaxis_EventRates_Elep_erec_afteradaptive_chain_1_job_0test_diagnostics.root",
        "/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/Onaxis_EventRates_Elep_erec_afteradaptive_chain_1_job_0test_diagnostics.root",
         "/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/Onaxis_EventRates_Elep_erec_afteradaptive_chain_1_job_0test_diagnostics.root",
         "/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/Onaxis_EventRates_Elep_erec_afteradaptive_chain_1_job_0test_diagnostics.root"
    };*/


std::vector<std::string> files_A = {//"/scratch/abipeake/Off_axis_chains_new/Septermber_q0q3_10percentonaxis/Septermber_q0q3_10percentonaxis_chain_0_job_0_diagnostics.root",
    "/scratch/abipeake/Off_axis_chains_new/seotember_offaxislongq0q3chain/seotember_offaxislongq0q3chain_chain_0_job_0_diagnostics.root",
   "/scratch/abipeake/Off_axis_chains_new/September_q0q3_75percentonaxis_part2/September_q0q3_75percentonaxis_part2_chain_0_job_0_diagnostics.root"};
std::vector<std::string> files_B = {//"/scratch/abipeake/Off_axis_chains_new/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive_chain_0_job_0copy_diagnostics.root",
"/scratch/abipeake/Off_axis_chains_new/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive_chain_0_job_0copy_diagnostics.root",
"/scratch/abipeake/Off_axis_chains_new/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive/LONGERWALLTIME2_q0q3_pTpzEnu_afteradaptive_chain_0_job_0copy_diagnostics.root"};


int main(){
// Open all files

    std::vector<TFile*> fA, fB;
    for (auto &fname : files_A) fA.push_back(TFile::Open(fname.c_str(), "READ"));
    for (auto &fname : files_B) fB.push_back(TFile::Open(fname.c_str(), "READ"));

    std::string pdfOut = "runplancomparisons_eleperec.pdf";
    TCanvas* c = new TCanvas("c","c",1200,2400);
    c->Print((pdfOut + "[").c_str());

    // Pattern for slice histograms
    std::regex eventPattern("^h_event_slice_[0-9]+$");

    // Labels for each set
std::vector<std::string> ratioLabels = {
    "10% on axis",
    //"25% on axis",
    "50% on axis",
    //"75% on axis",
    //"75% on axis"
};

// Colors for each set
const int colors[5] ={kRed, kOrange-3,kGreen+1, kCyan+2,kViolet+1};

// Loop over slices using the first file as reference
TIter nextkey(fA[0]->GetListOfKeys());
TKey* key;

while ((key = (TKey*)nextkey())) {
    std::string name = key->GetName();
    name = name.substr(0, name.find(";"));

    // Only process event slice histograms
    if (!std::regex_match(name, std::regex("^h_event_slice_[0-9]+$"))) continue;

    // Vectors for histograms
    std::vector<TH1D*> h_event_A_vec, h_event_B_vec;
    std::vector<TH1D*> h_post_A_vec, h_post_B_vec;

    bool missing = false;
    for (size_t i = 0; i < files_A.size(); ++i) {
        std::string postName = std::regex_replace(name, std::regex("h_event"), "h_post");

        TH1D* h_eA = (TH1D*) fA[i]->Get(name.c_str());
        TH1D* h_eB = (TH1D*) fB[i]->Get(name.c_str());
        TH1D* h_pA = (TH1D*) fA[i]->Get(postName.c_str());
        TH1D* h_pB = (TH1D*) fB[i]->Get(postName.c_str());

        if (!h_eA || !h_eB || !h_pA || !h_pB) {
            std::cerr << "Skipping slice " << name << " in set " << i
                      << " (missing histograms in " << files_A[i] << " or " << files_B[i] << ")\n";
            missing = true;
            break;
        }

        h_event_A_vec.push_back(h_eA);
        h_event_B_vec.push_back(h_eB);
        h_post_A_vec.push_back(h_pA);
        h_post_B_vec.push_back(h_pB);
    }

    if (missing) continue;

    std::string sliceLabel = h_event_A_vec[0]->GetTitle();
    double xmin = h_event_A_vec[0]->GetBinLowEdge(1);
    double xmax = h_event_A_vec[0]->GetBinLowEdge(h_event_A_vec[0]->GetNbinsX()) + 
                  h_event_A_vec[0]->GetBinWidth(h_event_A_vec[0]->GetNbinsX());

    c->Clear();
    c->Divide(1,3);

    const size_t nSets = h_event_A_vec.size(); // should be 3

    // -----------------------
    // Pad 1: Event Rates
    // -----------------------
    c->cd(1);
    
    TLegend* leg1 = new TLegend(0.55,0.7,0.88,0.88);

    for (size_t i = 0; i < nSets; ++i) {
        TH1D* hA = (TH1D*) h_event_A_vec[i]->Clone(Form("hA_%zu", i));
        hA->SetLineColor(colors[i]);
        hA->SetLineWidth(2);
        hA->SetMarkerStyle(0);

        if(i==0) hA->Draw("HIST E");
        else hA->Draw("HIST E SAME");

        leg1->AddEntry(hA, ratioLabels[i].c_str(), "l");
    }
    leg1->Draw();

        // Add labels directly on pad
        TLatex latex;
        latex.SetNDC(); latex.SetTextSize(0.03); latex.SetTextColor(kBlack);
        for (size_t i = 0; i < nSets; ++i)
            latex.DrawLatex(0.15, 0.85 - 0.05*i, ratioLabels[i].c_str());

    // -----------------------
    // Pad 2: Posterior uncertainties
    // -----------------------
    
    c->cd(2);
    TLegend* leg2 = new TLegend(0.15,0.7,0.45,0.88);

    for (size_t i = 0; i < nSets; ++i) {
        TH1D* hA = (TH1D*) h_post_A_vec[i]->Clone(Form("hPA_%zu", i));
        hA->SetLineColor(colors[i]);
        hA->SetLineWidth(2);
        hA->SetMarkerStyle(0);

        if(i==0) hA->Draw("HIST E");
        else hA->Draw("HIST E SAME");

        leg2->AddEntry(hA, ratioLabels[i].c_str(), "l");
    }
    leg2->Draw();

    // -----------------------
    // Pad 3: Posterior uncertainty ratios
    // -----------------------
    c->cd(3);
    double ymin = 1e6;
    double ymax = -1e6;

    for(size_t i=0; i<nSets; ++i){
        for(int b=1; b<=h_post_A_vec[i]->GetNbinsX(); ++b){
            double errA = h_post_A_vec[i]->GetBinError(b);
            double errB = h_post_B_vec[i]->GetBinError(b);
            if(errB==0) continue; // avoid division by zero
            double ratio = errA/errB;
            if(ratio < ymin) ymin = ratio;
            if(ratio > ymax) ymax = ratio;
        }
    }

// Optional: add some padding
double ypad = 0.1*(ymax - ymin);
ymin -= ypad;
ymax += ypad;

    TLegend* leg3 = new TLegend(0.15,0.7,0.45,0.88);

    // Draw horizontal dotted line at y=1
    TLine* line1 = new TLine(xmin, 1.0, xmax, 1.0);
    line1->SetLineColor(kBlack);
    line1->SetLineStyle(2); // dotted
    line1->SetLineWidth(2);
    line1->Draw();  // this draws first

    // Draw each set as smooth line
    for(size_t i=0; i<nSets; ++i){
        std::vector<double> x_vals, y_vals;
        for(int b=1;b<=h_post_A_vec[i]->GetNbinsX();++b){
            double binCenter = h_post_A_vec[i]->GetBinCenter(b);
            double errA = h_post_A_vec[i]->GetBinError(b);
            double errB = h_post_B_vec[i]->GetBinError(b);
            double ratio = (errB>0)? errA/errB : 0.0;
            x_vals.push_back(binCenter);
            y_vals.push_back(ratio);
        }

        TGraph* g = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());
        g->SetLineColor(colors[i]);
        g->SetLineWidth(2);
        g->SetMarkerStyle(0);
        g->SetTitle("Posterior Uncertainty Ratio");

        if(i==0){
            g->Draw("AL");      // Draw axes + line
            g->GetYaxis()->SetRangeUser(ymin, ymax); // scale y-axis to fit
        } else g->Draw("L SAME"); // subsequent lines
        leg3->AddEntry(g, ratioLabels[i].c_str(), "l");
    }

    leg3->Draw();

    // Slice label
    c->cd();
    latex.SetTextSize(0.04); latex.SetTextAlign(22);
    latex.DrawLatex(0.5, 0.96, sliceLabel.c_str());

    c->Print(pdfOut.c_str());
}

    c->Print((pdfOut + "]").c_str());

    for(auto f: fA) f->Close();
    for(auto f: fB) f->Close();
    delete c;

    return 0;
}