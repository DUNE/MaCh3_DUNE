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


#include <string>

#include <ctime>
#include <sstream>

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
void extract_2D_bins_from_yaml(const std::string& yaml_file, std::vector<BinDef>& binDefs, std::vector<double>& q0_edges, std::vector<double>& q3_edges, std::string xsecvar1, std::string xsecvar2) {
  std::set<double> q0_set;
  std::set<double> q3_set;
  binDefs.clear();

  try {
  YAML::Node config = YAML::LoadFile(yaml_file);
  const auto& systs = config["Systematics"];
  int index = 0;

  for (const auto& systematic : systs) {
  const auto& sys = systematic["Systematic"];
  const auto& cuts = sys["KinematicCuts"];

  double q0_min = -1, q0_max = -1, q3_min = -1, q3_max = -1;

  for (const auto& cut : cuts) {
  if (cut[xsecvar1.c_str()]) {
  const auto& q0 = cut[xsecvar1.c_str()];
  if (q0.IsSequence() && q0.size() == 2) {
  q0_min = q0[0].as<double>();
  q0_max = q0[1].as<double>();
  q0_set.insert(q0_min);
  q0_set.insert(q0_max);
  }
  }
  if (cut[xsecvar2.c_str()]) {
  const auto& q3 = cut[xsecvar2.c_str()];
  if (q3.IsSequence() && q3.size() == 2) {
  q3_min = q3[0].as<double>();
  q3_max = q3[1].as<double>();
  q3_set.insert(q3_min);
  q3_set.insert(q3_max);
  }
  }
  }
  if (std::isfinite(q0_min) ) {
      binDefs.push_back({index++, q0_min, q0_max});
  }
  if (std::isfinite(q3_min) ) {
      binDefs.push_back({index++, q3_min, q3_max});
  }
  }
  q0_edges.assign(q0_set.begin(), q0_set.end());
  q3_edges.assign(q3_set.begin(), q3_set.end());

  std::cout << "Parsed " << binDefs.size() << " bins from YAML\n";
  std::cout << "   q0 bins: " << q0_edges.size() - 1 << "\n";
  std::cout << "   q3 bins: " << q3_edges.size() - 1 << "\n";
  } catch (const std::exception& e) {
  std::cerr << " YAML parsing failed: " << e.what() << "\n";
  }
}

void extract_1D_bins_from_yaml(const std::string& yaml_file,
  std::vector<BinDef>& binDefs,
  std::vector<double>& q0_edges) {
  std::set<double> q0_set;
  binDefs.clear();

  try {
  YAML::Node config = YAML::LoadFile(yaml_file);
  const auto& systs = config["Systematics"];
  int index = 0;

  for (const auto& systematic : systs) {
  const auto& sys = systematic["Systematic"];
  const auto& cuts = sys["KinematicCuts"];

  double q0_min = -1, q0_max = -1, q3_min = -1, q3_max = -1;

  for (const auto& cut : cuts) {
  if (cut["RecoNeutrinoEnergy"]) {
  const auto& q0 = cut["RecoNeutrinoEnergy"];
  if (q0.IsSequence() && q0.size() == 2) {
  q0_min = q0[0].as<double>();
  q0_max = q0[1].as<double>();
  q0_set.insert(q0_min);
  q0_set.insert(q0_max);
  }
  }
  }
  if (std::isfinite(q0_min) ) {
      binDefs.push_back({index++, q0_min, q0_max});
  }
  }
  q0_edges.assign(q0_set.begin(), q0_set.end());

  std::cout << "Parsed " << binDefs.size() << " bins from YAML\n";
  std::cout << "   q0 bins: " << q0_edges.size() - 1 << "\n";
  } catch (const std::exception& e) {
  std::cerr << " YAML parsing failed: " << e.what() << "\n";
  }
}
*/

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

            for (const auto& cut : cuts) {

                std::cout << "Checking cut: ";
                if (cut[xsecvar1]) std::cout << xsecvar1 << " exists, ";
                else std::cout << xsecvar1 << " missing, ";
                if (cut[xsecvar2]) std::cout << xsecvar2 << " exists\n";
                else std::cout << xsecvar2 << " missing\n";

                if (cut[xsecvar1] && cut[xsecvar2]) {
                    const auto& q0 = cut[xsecvar1];
                    const auto& q3 = cut[xsecvar2];

                    if (q0.IsSequence() && q0.size() == 2 &&
                        q3.IsSequence() && q3.size() == 2) {

                        double q0_min = q0[0].as<double>();
                        double q0_max = q0[1].as<double>();
                        double q3_min = q3[0].as<double>();
                        double q3_max = q3[1].as<double>();

                        q0_set.insert(q0_min);
                        q0_set.insert(q0_max);
                        q3_set.insert(q3_min);
                        q3_set.insert(q3_max);

                        result.binDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
                    }
                }
            }
        }

        result.q0_edges.assign(q0_set.begin(), q0_set.end());
        result.q3_edges.assign(q3_set.begin(), q3_set.end());

        std::cout << "Parsed " << result.binDefs.size() << " 2D bins\n";
        std::cout << "   q0 edges: " << result.q0_edges.size() - 1 << "\n";
        std::cout << "   q3 edges: " << result.q3_edges.size() - 1 << "\n";

    } catch (const std::exception& e) {
        std::cerr << "YAML parsing failed: " << e.what() << "\n";
    }

    return result;
}
*/

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


void FixLowStatParams(covarianceXsec* xsec, TH2D* h_q0q3, double total_events, double frac_threshold,
                      std::vector<std::string>& fixed_names_out, const std::vector<BinDef>& binDefs) {
    int nbins_q0 = h_q0q3->GetNbinsY();
    int nbins_q3 = h_q0q3->GetNbinsX();
    double threshold = frac_threshold * total_events;

    std::cout << "Fixing low-stat parameters...\n";
    std::cout << " - frac_threshold (YAML) = " << frac_threshold << std::endl;
    std::cout << " - total events in histogram = " << total_events << std::endl;
    std::cout << " - Fixing threshold = " << threshold << " events\n";

    for (const auto& bin : binDefs) {
        double q0 = 0.5 * (bin.q0_min + bin.q0_max);
        double q3 = 0.5 * (bin.q3_min + bin.q3_max);

        int bin_q0 = h_q0q3->GetYaxis()->FindBin(q0);
        int bin_q3 = h_q0q3->GetXaxis()->FindBin(q3);

        if (bin_q0 < 1 || bin_q0 > nbins_q0 || bin_q3 < 1 || bin_q3 > nbins_q3) {
            std::cerr << "Warning: bin center out of range (q0=" << q0 << ", q3=" << q3 << ")\n";
            continue;
        }

        double val = h_q0q3->GetBinContent(bin_q3, bin_q0);
        std::cout << " - Bin idx=" << bin.index << ", content=" << val << " at (q0=" << q0 << ", q3=" << q3 << ")\n";

        if (val < threshold) {
            xsec->setSingleParameter(bin.index, 0.0);
            xsec->toggleFixParameter(bin.index); // If needed, guard with IsParameterFixed()
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

  std::string xsec_param_yaml = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_10.0GeV_fullgrid_smallerbins.yaml";
  auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
  auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
  auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();
  auto fixing_threshold = FitManager->raw()["General"]["Systematics"]["Parameter_fixing_threshold"].as<double>();

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

   
  //////////////////MaCh3 stuff
  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;
  std::vector<samplePDFFDBase*> DUNEPdfs;
 
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);
  std::unique_ptr<mcmc> MaCh3Fitter = std::make_unique<mcmc>(FitManager);
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

  TH2D* h_q0q3_combined = nullptr;
  std::cout << "DUNEPdfs size: " << DUNEPdfs.size() << std::endl;
 TFile* f = new TFile("debug_hist.root", "RECREATE");
  for (auto& pdf : DUNEPdfs) {
    auto* dunePdf = dynamic_cast<samplePDFDUNEBeamFD*>(pdf);
    if (!dunePdf) continue;
    std::vector<KinematicCut> SelectionVector;

    std::cout << "Calling get2DVarHist with vars: " << xsec_var1 << ", " << xsec_var2 << std::endl;
    std::cout << "AxisX bins: " << axisX->GetNbins()
          << ", AxisY bins: " << axisY->GetNbins() << std::endl;

    


    TH2* h = dunePdf->get2DVarHist(xsec_var1, xsec_var2, SelectionVector , /*WeightStyle=*/0, axisX, axisY);
    if (!h) {
        std::cerr << "Warning: get2DVarHist returned null.\n";
        continue;
    }
    std::cout << "Hist nbinsX = " << h->GetNbinsX() << ", axisX = " << axisX->GetNbins() << std::endl;

    std::cout << "Hist get2DVarHist integral: " << h->Integral() << std::endl;
    h->Write();
    if (!h_q0q3_combined) {
        h_q0q3_combined = (TH2D*)h->Clone("h_q0q3_combined");
        std::cout << "[DEBUG] Cloned histogram integral: " << h_q0q3_combined->Integral() << std::endl;

        h_q0q3_combined->SetDirectory(nullptr); // optional: detach from file
    } else {
        h_q0q3_combined->Add(h);
    }
}
f->Close();

  //  Print final integral
  if (h_q0q3_combined) {
      double h_q0q3_combined_int = h_q0q3_combined->Integral();
      std::cout << "Integral of h_q0q3_combined = " << h_q0q3_combined_int << std::endl;
  } else {
      std::cerr << "ERROR: No histograms were combined.\n";
  }

  std::vector<std::string> fixed_names;

  if (!StartFromPreviousChain) xsec->throwParameters();

  std::cout << "Fixed parameters due to low statistics:\n";
  for (const auto& name : fixed_names) {
      std::cout << "  - " << name << std::endl;
  }
  //FixLowStatParams(xsec, h_q0q3, total, 0.0001, fixed_names, binDefs);
  std::cout << "[DEBUG] Final combined hist integral: " << h_q0q3_combined->Integral() << std::endl;
for (int ix = 1; ix <= h_q0q3_combined->GetNbinsX(); ++ix) {
  for (int iy = 1; iy <= h_q0q3_combined->GetNbinsY(); ++iy) {
    double content = h_q0q3_combined->GetBinContent(ix, iy);
    if (content > 0)
      std::cout << "Bin (" << ix << ", " << iy << ") = " << content << std::endl;
  }
}

 std::vector<std::string> fixed_param_names;
 FixLowStatParams(xsec, h_q0q3_combined, h_q0q3_combined->Integral(), fixing_threshold, fixed_param_names, binning.binDefs);


  std::cout << "Fixed " << fixed_names.size() << " low-stat parameters." << std::endl;

  // === Record frozen parameters ===
  std::ofstream fixedOut("fixed_parameters.txt");
  for (const auto& name : fixed_names) fixedOut << name << "\n";
  fixedOut.close();


  // Re-freeze after throw
  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    if (std::find(fixed_names.begin(), fixed_names.end(), xsec->GetParFancyName(i)) != fixed_names.end()) {
      xsec->setSingleParameter(i, 1.0);
      if (!xsec->isParameterFixed(i)) xsec->toggleFixParameter(i);
    }
  }

  if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false))
    MaCh3Fitter->addSystObj(xsec);

  for (auto& pdf : DUNEPdfs) MaCh3Fitter->addSamplePDF(pdf);
  if (StartFromPreviousChain)
    MaCh3Fitter->StartFromPreviousFit(FitManager->raw()["General"]["PosFileName"].as<std::string>());

  
  // === Prepare Histograms ===
  TH1D* h_event_rate = (TH1D*)h_q0q3_combined->Clone("h_event_rate");
  TH1D* h_param_init  = (TH1D*)h_q0q3_combined->Clone("h_param_init");
  TH1D* h_param_frozen = (TH1D*)h_q0q3_combined->Clone("h_param_frozen");
  TH1D* h_param_mean   = (TH1D*)h_q0q3_combined->Clone("h_param_mean");
  TH1D* h_param_stddev = (TH1D*)h_q0q3_combined->Clone("h_param_stddev");

  h_param_init->Reset("ICES"); //remove the bin contents and stats error from histograms so I can refill them
  h_param_frozen->Reset("ICES");
  h_param_mean->Reset("ICES");
  h_param_stddev->Reset("ICES");

  h_event_rate->SetDirectory(nullptr);
  h_param_init->SetDirectory(nullptr);
  h_param_frozen->SetDirectory(nullptr);
  h_param_mean->SetDirectory(nullptr);
  h_param_stddev->SetDirectory(nullptr);

  std::cout << "[DEBUG] h_q0q3->IsOnHeap(): " << h_q0q3_combined->IsOnHeap() << std::endl;
  std::cout << "[DEBUG] gDirectory: " << gDirectory->GetName() << std::endl;
  


  MaCh3Fitter->runMCMC();


  /*
  std::cout << "[DEBUG] h_q0q3->IsOnHeap(): " << h_q0q3_combined->IsOnHeap() << std::endl;
  std::cout << "[DEBUG] gDirectory: " << gDirectory->GetName() << std::endl;

  std::cout << "done run MCMC" << std::endl;
  // === Load Posterior ===
  TFile* f_chain = new TFile(OutputFileName.c_str());
  TTree* post = (TTree*)f_chain->Get("posteriors");
  if (!post) {
    std::cerr << "Could not find 'posteriors' TTree in " << OutputFileName << std::endl;
    return 1;
  }
  if (!h_q0q3_combined) {
    std::cerr << "h_q0q3 is null before cloning!" << std::endl;
    return 1;
  }
  if (q0_edges.empty()) {
    std::cerr << " Error: q0 or q3 edges are empty. Histogram cannot be created." << std::endl;
    return 1;
  }
  //std::cout << "[DEBUG] About to clone h_q0q3: " << h_q0q3 << ", title: " << h_q0q3->GetTitle() << std::endl;


  // === Fill init/frozen state ===
  for (const auto& bin : binDefs) {
    double q0 = 0.5 * (bin.q0_min + bin.q0_max);
    
    int bin_q0 = h_q0q3_combined ->GetXaxis()->FindBin(q0);
    

    if (bin_q0 <= 0 || bin_q0 > h_q0q3_combined ->GetNbinsX()) continue;

    h_param_init->SetBinContent( bin_q0, 1.0);
    double frozen_val = xsec->isParameterFixed(bin.index) ? 0.0 : 1.0;
    h_param_frozen->SetBinContent( bin_q0, frozen_val);
  }

  // === Posterior mean/stddev ===
  int nParams = xsec->GetNumParams();
  std::vector<double> xsec_vals(nParams, 0.0);
  std::vector<double> sum(nParams, 0.0), sq_sum(nParams, 0.0);

  for (const auto& bin : binDefs) {
    post->SetBranchAddress(("xsec_" + std::to_string(bin.index)).c_str(), &xsec_vals[bin.index]);
  }

  int nEntries = post->GetEntries();
  for (int entry = 0; entry < nEntries; ++entry) {
    post->GetEntry(entry);
    for (const auto& bin : binDefs) {
      double val = xsec_vals[bin.index];
      sum[bin.index] += val;
      sq_sum[bin.index] += val * val;
    }
  }
  
  for (const auto& bin : binDefs) {
    double q0 = 0.5 * (bin.q0_min + bin.q0_max);
    int bin_q0 = h_q0q3_combined->GetXaxis()->FindBin(q0);

    double mean = sum[bin.index] / nEntries;
    double var  = std::max(0.0, (sq_sum[bin.index] / nEntries) - mean * mean);
    double stddev = std::sqrt(var);

    h_param_mean->SetBinContent(bin_q0, mean);
    h_param_stddev->SetBinContent(bin_q0, stddev);

    // Write to file
    outfile << bin.index << "\t" << q0 << "\t" << mean << "\t" << stddev << "\n";
  }

  // Close file
  outfile.close();

  // === Output diagnostics ===
  std::cout << "[CHECK] h_param_frozen integral: " << h_param_frozen->Integral() << std::endl;
  std::cout << "[CHECK] h_param_init integral: " << h_param_init->Integral() << std::endl;

  // === Save plots ===
  gStyle->SetPalette(kRainBow);
  TCanvas* c = new TCanvas("c", "", 900, 700);
  std::string pdfOutName = AddTimestampToFilename("param_summary_q0q3_someparamsfrozen1D.pdf");
  c->Print((pdfOutName + "[").c_str());  // Open multi-page PDF

  c->Clear(); h_event_rate->Draw("COLZ");     c->Print(pdfOutName.c_str());
  c->Clear(); h_param_init->Draw("COLZ");     c->Print(pdfOutName.c_str());
  c->Clear(); h_param_frozen->Draw("COLZ");   c->Print(pdfOutName.c_str());
  c->Clear(); h_param_mean->Draw("COLZ");     c->Print(pdfOutName.c_str());
  c->Clear(); h_param_stddev->Draw("COLZ");   c->Print(pdfOutName.c_str());

  c->Print((pdfOutName + "]").c_str());  // Close multi-page PDF


  //TFile* fOut = new TFile("param_summary_q0q3_someparamsfrozen.root", "RECREATE");
  std::string rootOutName = AddTimestampToFilename("param_summary_q0q3_someparamsfrozen1D.root");
  TFile* fOut = new TFile(rootOutName.c_str(), "RECREATE");

  h_q0q3_combined ->Write("h_q0q3_input");

  h_event_rate->Write();
  h_param_init->Write();
  h_param_frozen->Write();
  h_param_mean->Write();
  h_param_stddev->Write();
  fOut->Close();

  h_event_rate->Reset();
  h_param_init->Reset();
  h_param_frozen->Reset();
  h_param_mean->Reset();
  h_param_stddev->Reset();
  h_q0q3_combined ->Reset();


  std::cout << "\n✅ Saved:\n - param_summary_q0q3_someparamsfrozen.pdf\n - param_summary_q0q3_someparamsfrozen.root\n" << std::endl;
  std::cout << "Output ROOT file: " << rootOutName << std::endl;
  std::cout << "Output PDF file: " << pdfOutName << std::endl;
  
  */
}
