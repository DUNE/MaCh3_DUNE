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

std::string title;

// Forward declaration of your PDF type
class samplePDFFDBase;
struct PosteriorSample {
    std::vector<double> xsec_draw; // vector of xsec parameters for this posterior sample
};

// Example histogram axis struct
struct HistAxis {
    int bins;
    double min;
    double max;
    HistAxis(int b, double mn, double mx) : bins(b), min(mn), max(mx) {}
};


std::string AddTimestampToFilename(const std::string& baseName) {
    time_t now = time(0);
    tm* ltm = localtime(&now);
    std::ostringstream oss;
    oss << baseName.substr(0, baseName.find_last_of(".")) << "_"
        << 1900 + ltm->tm_year
        << 1 + ltm->tm_mon
        << ltm->tm_mday << "_"
        << ltm->tm_hour
        << ltm->tm_min
        << ltm->tm_sec << ".pdf";
    return oss.str();
}

BinningResult extract_2D_bins_from_yaml(const std::string& yaml_file, const std::string& xsecvar1, const std::string& xsecvar2) {
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
                if (cut[xsecvar1]) {
                    const auto& q0 = cut[xsecvar1];
                    if (q0.IsSequence() && q0.size() == 2) {
                        q0_min = q0[0].as<double>();
                        q0_max = q0[1].as<double>();
                        has_q0 = true;
                    }
                }
                if (cut[xsecvar2]) {
                    const auto& q3 = cut[xsecvar2];
                    if (q3.IsSequence() && q3.size() == 2) {
                        q3_min = q3[0].as<double>();
                        q3_max = q3[1].as<double>();
                        has_q3 = true;
                    }
                }
            }

            if (has_q0 && has_q3) {
                q0_set.insert(q0_min);
                q0_set.insert(q0_max);
                q3_set.insert(q3_min);
                q3_set.insert(q3_max);

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

  std::string formatAxisTitle(const std::string& var) {
    if (var == "q0") return "q_{0}";
    if (var == "q3") return "q_{3}";
    if (var == "TrueNeutrinoEnergy") return "E^_{#nu} [GeV]";
    if (var == "ENuProxy_minus_Enutrue") return "E^{True}_{#nu} -E^{Proxy}_{#nu} [GeV]";
    if (var == "yRec") return "y_{Rec}";
    if (var == "RecoNeutrinoEnergy") return "E^{Reco.}_{#nu} [GeV]";
    if (var == "isRelativeEnubias") return "E^{Bias}_{#nu}/E^_{#nu}  [GeV]";
    // Add more mappings as needed

    // Fallback: just return var unchanged if no match found
    return var;
}



int main(int argc, char* argv[]) {

    // --- Setup ROOT canvas ---
    TCanvas* c_master = new TCanvas("c", "c", 1200, 800);
    std::string pdfOut = AddTimestampToFilename("mcmc_summary_plots.pdf");

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mcmc_output.root> <config.yaml>" << std::endl;
        return 1;
    }

    std::string mcmc_file = argv[1];
    std::string config_file = argv[2];
    
    // -------------------------
    // Load YAML config
    // -------------------------
    YAML::Node config = YAML::LoadFile(config_file);
    manager* FitManager = new manager(config_file);
    if (!FitManager) {
        std::cerr << "Error: Failed to create FitManager from YAML config: " << config_file << std::endl;
        return 1;
    }

    // -------------------------
    // Derive output filenames
    // -------------------------
    std::string outputRootFile = mcmc_file.substr(0, mcmc_file.find_last_of('.')) + "_eventrate.root";
    std::string outputPdfFile  = mcmc_file.substr(0, mcmc_file.find_last_of('.')) + "_eventrate.pdf";

    auto OutputFile = std::unique_ptr<TFile>(TFile::Open(outputRootFile.c_str(), "RECREATE"));
    if (!OutputFile || OutputFile->IsZombie()) {
        std::cerr << "Error creating output file " << outputRootFile << std::endl;
        return 1;
    }
    OutputFile->cd();

    // -------------------------
    // Setup binning from YAML
    // -------------------------
    auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
    auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
    auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();

    if (xsec_yaml.empty()) {
        std::cerr << "Error: xsec_yaml path is empty or missing in config.\n";
        return 1;
    }

    BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);
    std::vector<BinDef> binDefs = binning.binDefs;
    std::vector<double> q0_edges = binning.q0_edges;
    std::vector<double> q3_edges = binning.q3_edges;

    if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Error: Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }

    TAxis* axisX = new TAxis(q0_edges.size() - 1, q0_edges.data());
    TAxis* axisY = new TAxis(q3_edges.size() - 1, q3_edges.data());

    // -------------------------
    // MaCh3 setup
    // -------------------------
    covarianceXsec* xsec = nullptr;
    covarianceOsc* osc = nullptr;
    std::vector<samplePDFFDBase*> DUNEPdfs;
    MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);
    std::unique_ptr<mcmc> MaCh3Fitter = std::make_unique<mcmc>(FitManager);

    // -------------------------
    // Build combined 2D histogram
    // -------------------------
    TH2D* h_q0q3_combined = nullptr;
    for (auto& pdf : DUNEPdfs) {
        auto* dunePdf = dynamic_cast<samplePDFDUNEBeamFD*>(pdf);
        if (!dunePdf) continue;

        std::vector<KinematicCut> SelectionVector;
        TH2* h = dunePdf->get2DVarHist(xsec_var1, xsec_var2, SelectionVector, 0, axisX, axisY);
        if (!h) {
            std::cerr << "Warning: get2DVarHist returned null.\n";
            continue;
        }
        if (!h_q0q3_combined) {
            h_q0q3_combined = (TH2D*)h->Clone("h_q0q3_combined");
            h_q0q3_combined->SetDirectory(nullptr);
        } else {
            h_q0q3_combined->Add(h);
        }
    }

    if (!h_q0q3_combined) {
        std::cerr << "Error: No 2D histograms were produced.\n";
        return 1;
    }

    // Write combined histogram to output diagnostics file
    h_q0q3_combined->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_q0q3_combined->GetYaxis()->SetTitle(xsec_var2.c_str());
    h_q0q3_combined->Write("xsec_eventrate_histo");

    // --- Event rate histogram with adaptive label colors ---
    TH2D* h_eventrate_labels = (TH2D*)h_q0q3_combined->Clone("h_eventrate_labels");
    h_eventrate_labels->SetTitle("Event Rate with Bin Labels;q_{0} [GeV];q_{3} [GeV]");

    //TCanvas* c_eventrate_labels = new TCanvas("c_eventrate_labels", "", 900, 700);
    //c_eventrate_labels->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kPastel);
    h_eventrate_labels->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_eventrate_labels->GetYaxis()->SetTitle(xsec_var2.c_str());

    h_eventrate_labels->Draw("COLZ");

    // Dynamic text size
    int nBinsX = h_eventrate_labels->GetNbinsX();
    int nBinsY = h_eventrate_labels->GetNbinsY();
    double textSize = 0.5 * std::min(1.0/nBinsX, 1.0/nBinsY);
    if (textSize > 0.05) textSize = 0.05;

    TLatex latex;
    latex.SetTextAlign(22); // center
    latex.SetTextSize(textSize);

    for (const auto& bin : binDefs) {
        double q0_center = 0.5 * (bin.q0_min + bin.q0_max);
        double q3_center = 0.5 * (bin.q3_min + bin.q3_max);

        int binX = h_eventrate_labels->GetXaxis()->FindBin(q0_center);
        int binY = h_eventrate_labels->GetYaxis()->FindBin(q3_center);
        double binContent = h_eventrate_labels->GetBinContent(binX, binY);

        // Adapt label color depending on bin content
        int color = kBlack;
        double maxContent = h_eventrate_labels->GetMaximum();
        if (binContent / maxContent > 0.5) color = kWhite;

        latex.SetTextColor(color);
        latex.DrawLatex(q0_center, q3_center, std::to_string(bin.index).c_str());
    }


        // --- Open multipage PDF ---
    //c_master->Print((outputFileName + "[").c_str()); 
    c_master->Print((outputPdfFile + "[").c_str());  // open PDF


    // 1️⃣ Event rate with labels
    c_master->Clear();
    h_eventrate_labels->Draw("COLZ TEXT");
    //c_master->Print(outputFileName.c_str());
    c_master->Print(outputPdfFile.c_str());          // add a page

    h_eventrate_labels->Write(); // save to ROOT

    // 2️⃣ Event rate heatmap without labels
    c_master->Clear();
    h_q0q3_combined->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_q0q3_combined->GetYaxis()->SetTitle(xsec_var2.c_str());

    h_q0q3_combined->Draw("COLZ");
    c_master->Update();
    //c_master->Print(outputFileName.c_str());
    c_master->Print(outputPdfFile.c_str());          // add a page


    // --- Close multipage PDF ---
    //c_master->Print((outputFileName + "]").c_str());
    c_master->Print((outputPdfFile + "]").c_str());  // close PDF




    // -------------------------
    // Clean up
    // -------------------------
    delete h_q0q3_combined;
    delete axisX;
    delete axisY;
    delete FitManager;


    OutputFile->Write();
    OutputFile->Close();

    std::cout << "Finished diagnostics. Output written to: " << outputPdfFile << std::endl;
    return 0;
}
