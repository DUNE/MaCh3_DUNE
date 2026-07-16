#include "Samples/BinningHandler.h"
#include "Samples/SampleHandlerBase.h"
#include "Samples/SampleHandlerFD.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Manager/MaCh3Logger.h"
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerBase.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <set>
#include <sstream>
#include <ctime>
#include <random>
#include <fstream>

#include "yaml-cpp/yaml.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"

struct BinDef {
    int index;
    double q0_min, q0_max;
    double q3_min, q3_max;
};

struct BinningResult {
    std::vector<BinDef> binDefs;
    std::vector<double> q0_edges;
    std::vector<double> q3_edges;
};

std::string AddTimestampToROOTFilename(const std::string& baseName) {
    time_t now = time(0);
    tm* ltm = localtime(&now);
    std::ostringstream oss;
    oss << baseName.substr(0, baseName.find_last_of(".")) << "_"
        << 1900 + ltm->tm_year
        << 1 + ltm->tm_mon
        << ltm->tm_mday << "_"
        << ltm->tm_hour
        << ltm->tm_min
        << ltm->tm_sec << ".root";
    return oss.str();
}

void setRedWhiteBluePalette() {
    const Int_t nRGBs = 3;
    Double_t stops[nRGBs] = {0.00, 0.50, 1.00};
    Double_t red[nRGBs]   = {0.00, 1.00, 1.00};
    Double_t green[nRGBs] = {0.00, 1.00, 0.00};
    Double_t blue[nRGBs]  = {1.00, 1.00, 0.00};
    const Int_t nColors = 255;
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nColors);
    gStyle->SetNumberContours(nColors);
}
void setOrangePalette() {
    const Int_t nRGBs = 5;
    Double_t stops[nRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
    Double_t red[nRGBs]   = { 1.00, 1.00, 0.95, 0.85, 0.55 };
    Double_t green[nRGBs] = { 0.98, 0.75, 0.50, 0.25, 0.00 };
    Double_t blue[nRGBs]  = { 0.92, 0.35, 0.10, 0.00, 0.00 };
    const Int_t nColors = 255;
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nColors);
    gStyle->SetNumberContours(nColors);
}
// -------------------------
// Helper: get latest non-empty posteriors tree
// -------------------------
TTree* getLatestNonEmpty(TFile* f) {
    int cycle = 9999;
    while (cycle > 0) {
        TString name = Form("posteriors;%d", cycle);
        TTree* t = (TTree*) f->Get(name);
        if (t && t->GetEntries() > 0) return t;
        cycle--;
    }
    return nullptr;
}

BinningResult extract_2D_bins_from_yaml(const std::string& yaml_file,
                                        const std::string& xsecvar1,
                                        const std::string& xsecvar2) {
    BinningResult result;
    std::set<double> q0_set;
    std::set<double> q3_set;

    try {
        YAML::Node config = YAML::LoadFile(yaml_file);
        const auto& systs = config["Systematics"];
        int index = 0;

        for (const auto& systematic : systs) {
            const auto& sys = systematic["Systematic"];
            const auto& cuts = sys["KinematicCuts"];

            double q0_min = 0, q0_max = 0;
            double q3_min = 0, q3_max = 0;
            bool has_q0 = false, has_q3 = false;

            for (const auto& cut : cuts) {
                if (cut[xsecvar1] && cut[xsecvar1].IsSequence() && cut[xsecvar1].size() == 2) {
                    double q0_lo = cut[xsecvar1][0].as<double>();
                    double q0_hi = cut[xsecvar1][1].as<double>();
                    q0_set.insert(q0_lo);
                    q0_set.insert(q0_hi);
                    q0_min = q0_lo;
                    q0_max = q0_hi;
                    has_q0 = true;
                }
                if (cut[xsecvar2] && cut[xsecvar2].IsSequence() && cut[xsecvar2].size() == 2) {
                    double q3_lo = cut[xsecvar2][0].as<double>();
                    double q3_hi = cut[xsecvar2][1].as<double>();
                    q3_set.insert(q3_lo);
                    q3_set.insert(q3_hi);
                    q3_min = q3_lo;
                    q3_max = q3_hi;
                    has_q3 = true;
                }
            }

            if (has_q0 || has_q3) {
                result.binDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
            }
        }

        result.q0_edges.assign(q0_set.begin(), q0_set.end());
        result.q3_edges.assign(q3_set.begin(), q3_set.end());

        std::cout << "Parsed " << result.binDefs.size() << " 2D bins\n";
        std::cout << "   q0 edges: " << (result.q0_edges.empty() ? 0 : result.q0_edges.size() - 1) << "\n";
        std::cout << "   q3 edges: " << (result.q3_edges.empty() ? 0 : result.q3_edges.size() - 1) << "\n";

    } catch (const std::exception& e) {
        std::cerr << "YAML parsing failed: " << e.what() << "\n";
    }

    return result;
}


int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mcmc_output.root> <config.yaml>" << std::endl;
        return 1;
    }

    std::string mcmc_file   = argv[1];
    std::string config_file = argv[2];

    std::string rootOut = AddTimestampToROOTFilename("mcmc_diagnostics.root");
    std::string pdfOut  = mcmc_file.substr(0, mcmc_file.find_last_of('.')) + "_diagnostics.pdf";

    // --- Open MCMC file ---
    TFile* f_post = TFile::Open(mcmc_file.c_str(), "READ");
    if (!f_post || f_post->IsZombie()) {
        std::cerr << "Cannot open MCMC file: " << mcmc_file << std::endl;
        return 1;
    }

    TTree* post = getLatestNonEmpty(f_post);
    if (!post) {
        std::cerr << "No non-empty 'posteriors' tree found in file." << std::endl;
        return 1;
    }
    std::cout << "Using " << post->GetName()
              << " with " << post->GetEntries() << " entries." << std::endl;

    // --- Load YAML config ---
    auto FitManager = std::make_unique<Manager>(config_file);
    auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
    auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
    auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();

    if (xsec_yaml.empty()) {
        std::cerr << "Error: xsec_yaml path is empty or missing in config.\n";
        return 1;
    }

    // --- Parse binning ---
    BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);
    std::vector<BinDef>& binDefs  = binning.binDefs;
    std::vector<double>& q0_edges = binning.q0_edges;
    std::vector<double>& q3_edges = binning.q3_edges;

    if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }

    TAxis* axisX = new TAxis(q0_edges.size() - 1, q0_edges.data());
    TAxis* axisY = new TAxis(q3_edges.size() - 1, q3_edges.data());

    // --- Setup MaCh3 PDFs ---
    ParameterHandlerGeneric* xsec = nullptr;
    std::vector<SampleHandlerFD*> DUNEPdfs;
    MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

    // --- Build 2D event-rate histogram ---
    TH2D* h_xsecvar_eventratehist = nullptr;
    for (auto& pdf : DUNEPdfs) {
        std::vector<KinematicCut> SelectionVector;
        TH2* h = pdf->Get2DVarHist(0, xsec_var1, xsec_var2, SelectionVector, 0, axisX, axisY);
        if (!h) {
            std::cerr << "Warning: Get2DVarHist returned null.\n";
            continue;
        }
        if (!h_xsecvar_eventratehist) {
            h_xsecvar_eventratehist = (TH2D*) h->Clone("h_xsecvar_eventratehist");
            h_xsecvar_eventratehist->SetDirectory(nullptr);
        } else {
            h_xsecvar_eventratehist->Add(h);
        }
    }

    if (!h_xsecvar_eventratehist) {
        std::cerr << "Error: No 2D histograms were produced.\n";
        return 1;
    }

    h_xsecvar_eventratehist->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_xsecvar_eventratehist->GetYaxis()->SetTitle(xsec_var2.c_str());

    // // --- Setup posterior parameter reading ---
    // int nParams = 0;
    // TObjArray* branches = post->GetListOfBranches();
    // for (int i = 0; i < branches->GetEntries(); ++i) {
    //     std::string name = branches->At(i)->GetName();
    //     if (name.rfind("param_", 0) == 0) {
    //         std::string suffix = name.substr(6);
    //         try {
    //             size_t pos;
    //             int idx = std::stoi(suffix, &pos);
    //             if (pos == suffix.size())
    //                 nParams = std::max(nParams, idx + 1);
    //         } catch (const std::invalid_argument&) {
    //             std::cout << "[Info] Skipping non-numeric param branch: " << name << "\n";
    //         }
    //     }
    // }
    // std::cout << "[Info] Found " << nParams << " param_* parameters\n";

    // std::vector<double> xsec_vals(nParams, 0.0);
    // post->SetBranchStatus("*", 0);
    // for (int i = 0; i < nParams; ++i) {
    //     std::string bname = "param_" + std::to_string(i);
    //     post->SetBranchStatus(bname.c_str(), 1);
    //     post->SetBranchAddress(bname.c_str(), &xsec_vals[i]);
    // }

    // --- Setup posterior parameter reading ---
    // Only read the param_N branches that correspond to parameters defined
    // in templateparams.yaml (i.e. those present in binDefs).
    int nParams = 0;
    for (const auto& bin : binDefs) {
        nParams = std::max(nParams, bin.index + 1);
    }
    std::cout << "[Info] Restricting to " << binDefs.size()
              << " parameters from templateparams.yaml (param_0 .. param_"
              << nParams - 1 << ")\n";

    std::vector<double> xsec_vals(nParams, 0.0);
    post->SetBranchStatus("*", 0);
    for (const auto& bin : binDefs) {
        std::string bname = "param_" + std::to_string(bin.index);
        if (!post->GetBranch(bname.c_str())) {
            std::cerr << "[Warning] Branch " << bname
                      << " not found in tree, skipping.\n";
            continue;
        }
        post->SetBranchStatus(bname.c_str(), 1);
        post->SetBranchAddress(bname.c_str(), &xsec_vals[bin.index]);
    }

    // --- Compute posterior mean and stddev per bin ---
    setRedWhiteBluePalette();
    TH2D* h_mean = new TH2D("h_param_mean", "Posterior Mean",
                            q0_edges.size()-1, q0_edges.data(),
                            q3_edges.size()-1, q3_edges.data());
    TH2D* h_stddev = new TH2D("h_param_stddev", "Posterior StdDev",
                              q0_edges.size()-1, q0_edges.data(),
                              q3_edges.size()-1, q3_edges.data());
    h_mean->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_mean->GetYaxis()->SetTitle(xsec_var2.c_str());
    h_stddev->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_stddev->GetYaxis()->SetTitle(xsec_var2.c_str());

    std::vector<double> sum(nParams, 0.0);
    std::vector<double> sq_sum(nParams, 0.0);

    Long64_t nEntries = post->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        post->GetEntry(i);
        for (const auto& bin : binDefs) {
            if (bin.index >= nParams) continue;
            double val = xsec_vals[bin.index];
            sum[bin.index]    += val;
            sq_sum[bin.index] += val * val;
        }
    }

    for (const auto& bin : binDefs) {
        if (bin.index >= nParams) continue;
        double mean   = sum[bin.index] / nEntries;
        double var    = std::max(0.0, (sq_sum[bin.index] / nEntries) - mean * mean);
        double stddev = std::sqrt(var);

        double q0_center = 0.5 * (bin.q0_min + bin.q0_max);
        double q3_center = 0.5 * (bin.q3_min + bin.q3_max);

        int bin_q0 = h_mean->GetXaxis()->FindBin(q0_center);
        int bin_q3 = h_mean->GetYaxis()->FindBin(q3_center);

        h_mean->SetBinContent(bin_q0, bin_q3, mean);
        h_stddev->SetBinContent(bin_q0, bin_q3, stddev);
    }

    // --- Write output ---
    auto OutputFile = std::unique_ptr<TFile>(TFile::Open(rootOut.c_str(), "RECREATE"));
    if (!OutputFile || OutputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << rootOut << std::endl;
        return 1;
    }
    OutputFile->cd();

    h_xsecvar_eventratehist->Write("xsec_eventrate_histo");
    h_mean->Write("xsec_param_mean");
    h_stddev->Write("xsec_param_stddev");

    // --- Save to PDF ---
    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    c->Print((pdfOut + "[").c_str());  // open PDF

    setRedWhiteBluePalette();
    h_xsecvar_eventratehist->Draw("COLZ");
    c->Print(pdfOut.c_str());
    c->Clear();

    h_mean->GetZaxis()->SetRangeUser(0.0, 4.0);
    h_mean->SetTitle("Posterior post-fit mean");
    h_mean->Draw("COLZ");
    c->Print(pdfOut.c_str());
    c->Clear();

    // With:
    setOrangePalette();                              // <-- switch palette
    h_stddev->GetZaxis()->SetRangeUser(0.0, 0.45);
    h_stddev->SetTitle("Posterior post-fit std dev");
    h_stddev->Draw("COLZ");
    c->Print(pdfOut.c_str());
    c->Clear();
    setRedWhiteBluePalette();                        // <-- restore for any later plots

    c->Print((pdfOut + "]").c_str());  // close PDF

    // --- Cleanup ---
    delete h_xsecvar_eventratehist;
    delete h_mean;
    delete h_stddev;
    delete axisX;
    delete axisY;
    delete c;
    f_post->Close();
    delete f_post;

    OutputFile->Write();
    OutputFile->Close();

    std::cout << "Saved diagnostics to: " << rootOut << std::endl;
    std::cout << "PDF written to: " << pdfOut << std::endl;
    return 0;
}
