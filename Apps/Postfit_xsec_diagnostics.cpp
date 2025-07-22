#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <ctime>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TROOT.h>

#include <yaml-cpp/yaml.h>

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

void Save1DSlices(TH2D* h2d, bool slice_along_q3, const std::string& label, TCanvas* c, const std::string& pdfOut) {
    int nSlices = slice_along_q3 ? h2d->GetNbinsY() : h2d->GetNbinsX();

    for (int i = 1; i <= nSlices; ++i) {
        std::string slice_name = Form("%s_slice_%d", h2d->GetName(), i);
        TH1D* h_slice = slice_along_q3 ? h2d->ProjectionX(slice_name.c_str(), i, i)
                                       : h2d->ProjectionY(slice_name.c_str(), i, i);

        // Set axis labels
        h_slice->SetTitle(Form("%s Slice %d", label.c_str(), i));
        h_slice->GetXaxis()->SetTitle(slice_along_q3 ? "E_{rec} - E_{#nu}" : "E_{#nu} [GeV]");
        h_slice->GetYaxis()->SetTitle("Posterior Value");

        h_slice->Draw("HIST");
        c->Print(pdfOut.c_str());
    }
}



int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mcmc_output.root> <bin_config.yaml>" << std::endl;
        return 1;
    }

    std::string mcmc_file = argv[1];
    std::string yaml_file = argv[2];

    TFile* fout = new TFile(rootOut.c_str(), "RECREATE");

    covarianceXsec* xsec = nullptr;
    covarianceOsc* osc = nullptr;
    std::vector<samplePDFFDBase*> DUNEPdfs;
    
    MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);

    // Load binning
    std::string xsec_param_yaml = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_10.0GeV_fullgrid_smallerbins.yaml";
    auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
    auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
    auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();

  /*
  std::ifstream test("/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_10.0GeV_fullgrid_smallerbins.yaml");
    if (!test.is_open()) {
    std::cerr << "File not found: /scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_10.0GeV_fullgrid_smallerbins.yaml" << std::endl;
    return 1;
  }*/
  //extract_2D_bins_from_yaml(xsec_yaml, binDefs, q0_edges, q3_edges, xsec_var1, xsec_var2);

  BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);

  TAxis* axisX = new TAxis(binning.q0_edges.size() - 1, binning.q0_edges.data());
  TAxis* axisY = new TAxis(binning.q3_edges.size() - 1, binning.q3_edges.data());
  
    // Lets draw the event rate distribution in terms of the xsec variables...


  for (auto& pdf : DUNEPdfs) {
    auto* dunePdf = dynamic_cast<samplePDFDUNEBeamFD*>(pdf);
    if (!dunePdf) continue;
    std::vector<KinematicCut> SelectionVector;

    TH2* h = dunePdf->get2DVarHist(xsec_var1, xsec_var2, SelectionVector , /*WeightStyle=*/0, axisX, axisY);
    if (!h) {
        std::cerr << "Warning: get2DVarHist returned null.\n";
        continue;
    }

    if (!h_q0q3_combined) {
        h_q0q3_combined = (TH2D*)h->Clone("h_q0q3_combined");
        h_q0q3_combined->SetDirectory(nullptr); // optional: detach from file
    } else {
        h_q0q3_combined->Add(h);
    }
}


    if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Error: Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }

    h->GetXaxis()->SetTitle(xsec_var1.c_str());
    h->GetYaxis()->SetTitle(xsec_var2.c_str());
    h->Write("xsec_eventrate_histo");
    fout->Close();

    // Load posterior
    TFile* f_post = new TFile(mcmc_file.c_str());
    TTree* post = (TTree*)f_post->Get("posteriors");
    if (!post) {
        std::cerr << "Error: Could not find 'posteriors' TTree in file " << mcmc_file << std::endl;
        return 1;
    }

    int nParams = 0;
    for (const auto& bin : binDefs)
        if (bin.index >= nParams)
            nParams = bin.index + 1;

    std::vector<double> xsec_vals(nParams, 0.0);
    std::vector<double> sum(nParams, 0.0), sq_sum(nParams, 0.0);

    for (const auto& bin : binDefs) {
        post->SetBranchAddress(("xsec_" + std::to_string(bin.index)).c_str(), &xsec_vals[bin.index]);
    }

    int nEntries = post->GetEntries();
    for (int i = 0; i < nEntries; ++i) {
        post->GetEntry(i);
        for (const auto& bin : binDefs) {
            double val = xsec_vals[bin.index];
            sum[bin.index] += val;
            sq_sum[bin.index] += val * val;
        }
    }


    TH2D* h_mean   = new TH2D("h_param_mean", "Posterior Mean;E_{rec} - E_{#nu} [GeV];True E_{#nu} [GeV]", q3_edges.size()-1, &q3_edges[0], q0_edges.size()-1, &q0_edges[0]);
    TH2D* h_stddev = new TH2D("h_param_stddev", "Posterior StdDev;E_{rec} - E_{#nu} [GeV];True E_{#nu} [GeV]", q3_edges.size()-1, &q3_edges[0], q0_edges.size()-1, &q0_edges[0]);

    for (const auto& bin : binDefs) {
        double mean = sum[bin.index] / nEntries;
        double var = std::max(0.0, (sq_sum[bin.index] / nEntries) - mean * mean);
        double stddev = std::sqrt(var);

        double q0 = 0.5 * (bin.q0_min + bin.q0_max);
        double q3 = 0.5 * (bin.q3_min + bin.q3_max);

        int bin_q0 = h_mean->GetYaxis()->FindBin(q0);
        int bin_q3 = h_mean->GetXaxis()->FindBin(q3);

        h_mean->SetBinContent(bin_q3, bin_q0, mean);
        h_stddev->SetBinContent(bin_q3, bin_q0, stddev);
    }

   
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetRightMargin(0.15);

    std::string pdfOut = AddTimestampToFilename("mcmc_summary_plots.pdf");
    c->Print((pdfOut + "[").c_str());

    c->Clear(); h_mean->Draw("COLZ"); c->Print(pdfOut.c_str());
    c->Clear(); h_stddev->Draw("COLZ"); c->Print(pdfOut.c_str());

    // Slice vertically: q3 fixed, q0 varying (i.e., slices along q3 axis)
    Save1DSlices(h_mean, true, "Posterior Mean", c, pdfOut);
    Save1DSlices(h_stddev, true, "Posterior StdDev", c, pdfOut);

    //Also slice horizontally (q0 fixed, q3 varying)
    Save1DSlices(h_mean, false, "Posterior Mean (q0 slices)", c, pdfOut);
    Save1DSlices(h_stddev, false, "Posterior StdDev (q0 slices)", c, pdfOut);


    c->Print((pdfOut + "]").c_str());

    std::cout << "Plots saved to: " << pdfOut << std::endl;

   
    std::string rootOut = AddTimestampToFilename("mcmc_summary_plots.root");
    rootOut.replace(rootOut.find(".pdf"), 4, ".root");
    
    h_mean->Write();
    h_stddev->Write();
    fout->Close();
    std::cout << "Histograms saved to: " << rootOut << std::endl;

    return 0;
}
