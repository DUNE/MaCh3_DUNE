#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <cmath> 

#include "mcmc/mcmc.h"
#include "samplePDFDUNE/MaCh3DUNEFactory.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>

#include <yaml-cpp/yaml.h>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem> 

#include <string>

#include <ctime>
#include <sstream>

namespace fs = std::filesystem; 

std::string AddTimestampToFilename(const std::string& baseName) { /////////////function to make sure I dont overwrite any of the files produced :)
    time_t now = time(0);
    tm *ltm = localtime(&now);
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
/*
// ADD THIS
struct BinDef {
  int index;
  double q0_min, q0_max;
  double q3_min, q3_max;
};*/

struct BinDef {
    int index;
    double q0_min;
    double q0_max;
    double q3_min;
    double q3_max;
};

struct BinningResult {
    std::vector<BinDef> binDefs;
    std::vector<double> q0_edges;
    std::vector<double> q3_edges;
};

/*
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
}*/

BinningResult extract_bins_from_yaml(const std::string& yaml_file,
                                     const std::string& xsecvar1,
                                     const std::string& xsecvar2) {
    BinningResult result;
    std::set<double> var1_edges;
    std::set<double> var2_edges;

    bool has_var1 = false;
    bool has_var2 = false;

    try {
        YAML::Node config = YAML::LoadFile(yaml_file);
        const auto& systs = config["Systematics"];
        int index = 0;

        for (const auto& systematic : systs) {
            const auto& sys = systematic["Systematic"];
            const auto& cuts = sys["KinematicCuts"];

            double v1_min = 0, v1_max = 0;
            double v2_min = 0, v2_max = 0;
            bool found_v1 = false, found_v2 = false;

            for (const auto& cut : cuts) {
                if (cut[xsecvar1]) {
                    const auto& v1 = cut[xsecvar1];
                    if (v1.IsSequence() && v1.size() == 2) {
                        v1_min = v1[0].as<double>();
                        v1_max = v1[1].as<double>();
                        found_v1 = true;
                        has_var1 = true;
                    }
                }
                if (cut[xsecvar2]) {
                    const auto& v2 = cut[xsecvar2];
                    if (v2.IsSequence() && v2.size() == 2) {
                        v2_min = v2[0].as<double>();
                        v2_max = v2[1].as<double>();
                        found_v2 = true;
                        has_var2 = true;
                    }
                }
            }

            if (found_v1 && found_v2) {
                // 2D bin
                var1_edges.insert(v1_min);
                var1_edges.insert(v1_max);
                var2_edges.insert(v2_min);
                var2_edges.insert(v2_max);
                result.binDefs.push_back({index++, v1_min, v1_max, v2_min, v2_max});
            } else if (found_v1 && !found_v2) {
                // 1D bin (var1 only)
                var1_edges.insert(v1_min);
                var1_edges.insert(v1_max);
                result.binDefs.push_back({index++, v1_min, v1_max, 0.0, 0.0});
            } else if (found_v2 && !found_v1) {
                // 1D bin (var2 only)
                var2_edges.insert(v2_min);
                var2_edges.insert(v2_max);
                result.binDefs.push_back({index++, v2_min, v2_max, 0.0, 0.0});
            }
        }

        // Assign edges based on what we found
        if (has_var1) result.q0_edges.assign(var1_edges.begin(), var1_edges.end());
        if (has_var2) result.q3_edges.assign(var2_edges.begin(), var2_edges.end());

        // Logging
        if (has_var1 && has_var2) {
            std::cout << "Parsed " << result.binDefs.size() << " 2D bins\n";
            std::cout << "   " << xsecvar1 << " edges: " 
                      << (result.q0_edges.empty() ? 0 : result.q0_edges.size() - 1) << "\n";
            std::cout << "   " << xsecvar2 << " edges: "
                      << (result.q3_edges.empty() ? 0 : result.q3_edges.size() - 1) << "\n";
        } else if (has_var1) {
            std::cout << "Parsed " << result.binDefs.size() << " 1D bins (variable: " 
                      << xsecvar1 << ")\n";
            std::cout << "   " << xsecvar1 << " edges: "
                      << (result.q0_edges.empty() ? 0 : result.q0_edges.size() - 1) << "\n";
        } else if (has_var2) {
            std::cout << "Parsed " << result.binDefs.size() << " 1D bins (variable: " 
                      << xsecvar2 << ")\n";
            std::cout << "   " << xsecvar2 << " edges: "
                      << (result.q3_edges.empty() ? 0 : result.q3_edges.size() - 1) << "\n";
        } else {
            std::cerr << "No valid binning variables (" << xsecvar1 << ", " << xsecvar2
                      << ") found in YAML.\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "YAML parsing failed: " << e.what() << "\n";
    }

    return result;
}



void FixLowStatParams(
    covarianceXsec* xsec,
    TH1* h_1d,             // can be TH1D or nullptr
    TH2D* h_2d,            // can be nullptr
    double total_events,
    double frac_threshold,
    std::vector<std::string>& fixed_names_out,
    const std::vector<BinDef>& binDefs
) {
    double threshold = frac_threshold * total_events;

    std::cout << "\n[FixLowStatParams] Starting low-stat parameter fixing...\n";
    std::cout << " - frac_threshold = " << frac_threshold << "\n";
    std::cout << " - total events = " << total_events << "\n";
    std::cout << " - event threshold = " << threshold << "\n";

    bool is2D = (h_2d != nullptr);
    bool is1D = (h_1d != nullptr);

    if (!is1D && !is2D) {
        std::cerr << "[ERROR] Both h_1d and h_2d are null — nothing to do.\n";
        return;
    }

    for (const auto& bin : binDefs) {
        double val = 0.0;
        int binx = -1, biny = -1;

        if (is2D) {
            double q0 = 0.5 * (bin.q0_min + bin.q0_max);
            double q3 = 0.5 * (bin.q3_min + bin.q3_max);

            binx = h_2d->GetXaxis()->FindBin(q0);
            biny = h_2d->GetYaxis()->FindBin(q3);
            val = h_2d->GetBinContent(binx, biny);

            std::cout << "Bin idx=" << bin.index
                      << " (q0=" << q0 << ", q3=" << q3
                      << ") → content=" << val << "\n";
        }
        else if (is1D) {
            // Pick the valid variable from the binDef
            double q = (bin.q0_max > bin.q0_min) ? 0.5 * (bin.q0_min + bin.q0_max)
                                                 : 0.5 * (bin.q3_min + bin.q3_max);
            binx = h_1d->GetXaxis()->FindBin(q);
            val = h_1d->GetBinContent(binx);

            std::cout << "Bin idx=" << bin.index
                      << " (center=" << q << ") → content=" << val << "\n";
        }

        // Out-of-range protection
        if (val < 0) val = 0.0;
        if (val < threshold) {
            std::cout << std::fixed << std::setprecision(2)
                      << "  → FIXING param " << bin.index
                      << " (val=" << val << " < threshold=" << threshold << ")\n";

            xsec->setSingleParameter(bin.index, 0.0);
            xsec->toggleFixParameter(bin.index);
            fixed_names_out.push_back(xsec->GetParFancyName(bin.index));
        }
    }

    std::cout << "Total fixed parameters: " << fixed_names_out.size() << "\n";
}




/// --- Main Fit and Diagnostics ---
int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: Fit config.yaml" << std::endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);


  //////// Open output file for mean and std
  std::ofstream outfile("param_stats.txt");
  outfile << "bin_index\tq0_center\tmean\tstddev\n";


  // === MCMC setup ===
  std::vector<BinDef> binDefs;
  std::vector<double> q0_edges;
  std::vector<double> q3_edges;
  
  manager* FitManager = new manager(argv[1]);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  
  if (!OutputFile || OutputFile->IsZombie()) {
    std::cerr << "Failed to open output file: " << OutputFileName << std::endl;
    return 1;
  }
  OutputFile->cd();

  
  auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
  auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
  auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();
  auto fixing_threshold = FitManager->raw()["General"]["Systematics"]["Parameter_fixing_threshold"].as<double>();

  
  //extract_2D_bins_from_yaml(xsec_yaml, binDefs, q0_edges, q3_edges, xsec_var1, xsec_var2);
  BinningResult binning = extract_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);

  //TAxis* axisX = new TAxis(binning.q0_edges.size() - 1, binning.q0_edges.data());
  //TAxis* axisY = new TAxis(binning.q3_edges.size() - 1, binning.q3_edges.data());
  // Determine dimensionality
bool is2D = (!binning.q0_edges.empty() && !binning.q3_edges.empty());
bool is1D_q0 = (!binning.q0_edges.empty() && binning.q3_edges.empty());
bool is1D_q3 = (!binning.q3_edges.empty() && binning.q0_edges.empty());

TAxis* axisX = nullptr;
TAxis* axisY = nullptr;

if (is2D) {
    axisX = new TAxis(binning.q0_edges.size() - 1, binning.q0_edges.data());
    axisY = new TAxis(binning.q3_edges.size() - 1, binning.q3_edges.data());
    std::cout << "Created 2D binning axes for " << xsec_var1 << " and " << xsec_var2 << "\n";
} 
else if (is1D_q0) {
    axisX = new TAxis(binning.q0_edges.size() - 1, binning.q0_edges.data());
    std::cout << "Created 1D binning axis for " << xsec_var1 << "\n";
}
else if (is1D_q3) {
    axisY = new TAxis(binning.q3_edges.size() - 1, binning.q3_edges.data());
    std::cout << "Created 1D binning axis for " << xsec_var2 << "\n";
}
else {
    std::cerr << "Error: No valid binning edges found in YAML for "
              << xsec_var1 << " or " << xsec_var2 << "\n";
    return 1;
}
   
  //////////////////MaCh3 stuff
  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;
  std::vector<samplePDFFDBase*> DUNEPdfs;
 
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);
  std::unique_ptr<mcmc> MaCh3Fitter = std::make_unique<mcmc>(FitManager);


  //stuff added from Liban

  std::string throwmatrixfilename = GetFromManager<std::string>(FitManager->raw()["General"]["ThrowMatrixFile"], "");
  std::string throwmatrixname = GetFromManager<std::string>(FitManager->raw()["General"]["ThrowMatrixName"], "");
  if (throwmatrixfilename == "") {
    MACH3LOG_INFO("No throw matrix file specified, will throw from covariance matrix.");
  }
  else {
    TFile *throwmatrixfile = new TFile(throwmatrixfilename.c_str());
    if (throwmatrixfile->IsZombie()) {
      MACH3LOG_ERROR("Couldn't find {}", throwmatrixfilename);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    MACH3LOG_INFO("Found throw matrix file {}.", throwmatrixfilename);
    TMatrixDSym *throwmatrix = throwmatrixfile->Get<TMatrixDSym>(throwmatrixname.c_str());
    if (!throwmatrix) {
      MACH3LOG_ERROR("Couldn't find throw matrix {} in file {}", throwmatrixname, throwmatrixfilename);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    // Reset individual step scales to 1.0
    xsec->resetIndivStepScale();
    // Keep global step scale flexible, but it should be 1
    double globalStepScale = FitManager->raw()["General"]["Systematics"]["XsecStepScale"].as<double>();
    if (globalStepScale != 1.0) {
      MACH3LOG_WARN("Global step scale is not 1.0, it is set to {}. This may cause issues when using adapted throw matrix.", globalStepScale);
    }
    xsec->setStepScale(globalStepScale);
    MACH3LOG_WARN("I have set all the individual step scales to 1.0 since we are using an external throw matrix");
    xsec->setThrowMatrix(throwmatrix);
    MACH3LOG_INFO("Set throw matrix from file {} with name {}",
                  throwmatrixfilename, throwmatrixname);
    // Print the throw matrix diagonals
    for (int i = 0; i < throwmatrix->GetNrows(); i++) {
      std::cout << (*throwmatrix)(i, i) << " ";
    }
    std::cout << std::endl;
  }

  //
  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);

  
  //Some place to store the histograms
  std::vector<TH1*> PredictionHistograms;
  std::vector<std::string> sample_names;

  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {
    
    std::string name = DUNEPdfs[sample_i]->GetTitle();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    
    osc->setParameters();
    DUNEPdfs[sample_i] -> reweight();
    if (DUNEPdfs[sample_i]->GetNDim() == 1){
      PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc")));
      DUNEPdfs[sample_i]->addData(static_cast<TH1D*>(PredictionHistograms[sample_i]));
    }
      
    else if (DUNEPdfs[sample_i]->GetNDim() == 2){
      PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i] -> get2DHist() -> Clone(NameTString+"_unosc")));
      DUNEPdfs[sample_i]->addData(static_cast<TH2D*>(PredictionHistograms[sample_i]));
    }
    else {
      MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting"); 
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
    
    
  }

  // === Create output file ===
TFile* f = new TFile("debug_hist.root", "RECREATE");
if (!f || f->IsZombie()) {
    std::cerr << "Failed to open debug_hist.root for writing.\n";
    return 1;
}


std::cout << "Detected binning mode: " 
          << (is2D ? "2D" : (is1D_q0 ? "1D (xsec_var1)" : "1D (xsec_var2)")) 
          << std::endl;

// Combined histograms (dimension-specific)
TH1* h_combined_1D = nullptr;
TH2D* h_combined_2D = nullptr;

std::cout << "DUNEPdfs size: " << DUNEPdfs.size() << std::endl;

for (auto& pdf : DUNEPdfs) {
    auto* dunePdf = dynamic_cast<samplePDFDUNEBeamFD*>(pdf);
    if (!dunePdf) continue;
    std::vector<KinematicCut> SelectionVector;

    // ---------------------
    // 2D binning case
    // ---------------------
      
    if (is2D) {
        TH2* h2 = dunePdf->get2DVarHist(xsec_var1, xsec_var2,
                                        SelectionVector, /*WeightStyle=*/0,
                                        axisX, axisY);
        if (!h2) {
            std::cerr << "[WARN] get2DVarHist returned null.\n";
            continue;
        }

        h2->GetXaxis()->SetTitle(xsec_var1.c_str());
        h2->GetYaxis()->SetTitle(xsec_var2.c_str());
        h2->Write("", TObject::kOverwrite);

        if (!h_combined_2D) {
            h_combined_2D = (TH2D*)h2->Clone("h_q0q3_combined");
            h_combined_2D->SetDirectory(nullptr);
        } else {
            h_combined_2D->Add(h2);
        }

        std::cout << "[DEBUG] Added 2D hist integral: " << h2->Integral() << std::endl;
    }

    // ---------------------
    // 1D binning case
    // ---------------------
    else if (is1D_q0 || is1D_q3) {
        std::string var = is1D_q0 ? xsec_var1 : xsec_var2;
        TAxis* axis = is1D_q0 ? axisX : axisY;

        TH1* h1 = dunePdf->get1DVarHist(var, SelectionVector, /*WeightStyle=*/0, axis);
        if (!h1) {
            std::cerr << "[WARN] get1DVarHist returned null.\n";
            continue;
        }

        h1->GetXaxis()->SetTitle(var.c_str());
        h1->Write("", TObject::kOverwrite);

        if (!h_combined_1D) {
            h_combined_1D = (TH1*)h1->Clone("h_combined_1D");
            h_combined_1D->SetDirectory(nullptr);
        } else {
            h_combined_1D->Add(h1);
        }

        std::cout << "[DEBUG] Added 1D hist integral: " << h1->Integral() << std::endl;
    }
}

f->Close();
std::cout << "Output file closed: debug_hist.root\n";

// ---------------------
// Final integral summary
// ---------------------
if (is2D && h_combined_2D) {
    double total = h_combined_2D->Integral();
    std::cout << "Final 2D combined histogram integral = " << total << std::endl;

    // Optional: print nonzero bins
    for (int ix = 1; ix <= h_combined_2D->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h_combined_2D->GetNbinsY(); ++iy) {
            double content = h_combined_2D->GetBinContent(ix, iy);
            if (content > 0)
                std::cout << "Bin (" << ix << ", " << iy << ") = " << content << std::endl;
        }
    }
}
else if (!is2D && h_combined_1D) {
    double total = h_combined_1D->Integral();
    std::cout << "Final 1D combined histogram integral = " << total << std::endl;

    // Optional: print nonzero bins
    for (int ix = 1; ix <= h_combined_1D->GetNbinsX(); ++ix) {
        double content = h_combined_1D->GetBinContent(ix);
        if (content > 0)
            std::cout << "Bin " << ix << " = " << content << std::endl;
    }
}
else {
    std::cerr << "ERROR: No histograms were combined.\n";
}

 std::vector<std::string> fixed_names;

if (is2D && h_combined_2D) {
    FixLowStatParams(xsec, nullptr, h_combined_2D, 
                     h_combined_2D->Integral(), fixing_threshold, 
                     fixed_names, binning.binDefs);
}
else if (!is2D && h_combined_1D) {
    FixLowStatParams(xsec, h_combined_1D, nullptr, 
                     h_combined_1D->Integral(), fixing_threshold, 
                     fixed_names, binning.binDefs);
}

  if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false))
    MaCh3Fitter->addSystObj(xsec);

  for (auto& pdf : DUNEPdfs) MaCh3Fitter->addSamplePDF(pdf);
  if (StartFromPreviousChain)
    MaCh3Fitter->StartFromPreviousFit(FitManager->raw()["General"]["PosFileName"].as<std::string>());

  
  
  MaCh3Fitter->runMCMC();

}
