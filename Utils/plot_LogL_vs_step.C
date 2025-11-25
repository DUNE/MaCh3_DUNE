#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TBranch.h>
#include <iostream>
#include <vector>

void plot_LogL_vs_step(const char* filename = "chain_output.root") {
    // Extract base name
    TString fullPath(filename);
    TString baseName = gSystem->BaseName(fullPath);
    baseName.ReplaceAll(".root", "");
    TString outPdf = baseName + "_LogL_and_traces.pdf";

    // Open ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Get tree
    TTree *tree = (TTree*)file->Get("posteriors");
    if (!tree) {
        std::cerr << "TTree 'posteriors' not found!" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Find first three parameters (skip step & LogL)
    std::vector<TString> paramNames;
    TObjArray *branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TString bname = branches->At(i)->GetName();
        if (bname == "step" || bname == "LogL") continue;
        paramNames.push_back(bname);
        if (paramNames.size() >= 3) break;
    }
    if (paramNames.size() < 3) {
        std::cerr << "Less than 3 parameters found in 'posteriors'!" << std::endl;
    }

    // Prepare TGraphs
    tree->Draw("LogL:step", "", "goff");
    TGraph *gLogL = new TGraph(tree->GetSelectedRows(),
                               tree->GetV2(), tree->GetV1());
    gLogL->SetTitle("LogL vs Step;Step;LogL");

    std::vector<TGraph*> gParams;
    for (auto &p : paramNames) {
        TString drawCmd = p + ":step";
        tree->Draw(drawCmd, "", "goff");
        gParams.push_back(new TGraph(tree->GetSelectedRows(),
                                     tree->GetV2(), tree->GetV1()));
        gParams.back()->SetTitle(Form("%s vs Step;Step;%s", p.Data(), p.Data()));
    }

    // Create multi-pad canvas
    int npads = 1 + gParams.size();
    TCanvas *c = new TCanvas("c", "LogL and Parameter Traces", 1000, 250*npads);
    c->Divide(1, npads);

    // Plot LogL
    c->cd(1);
    gLogL->Draw("AL");

    // Plot parameters
    for (size_t i = 0; i < gParams.size(); ++i) {
        c->cd(i+2);
        gParams[i]->Draw("AL");
    }

    // Save PDF
    c->SaveAs(outPdf);
    std::cout << "Saved plot as " << outPdf << std::endl;

    file->Close();
    delete file;
}
