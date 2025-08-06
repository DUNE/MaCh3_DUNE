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
    if (var == "TrueNeutrinoEnergy") return "E^{True}_{#nu} [GeV]";
    if (var == "ENuProxy_minus_Enutrue") return "E^{True}_{#nu} -E^{Proxy}_{#nu} [GeV]";
    if (var == "yRec") return "y_{Rec}";
    if (var == "RecoNeutrinoEnergy") return "E^{Reco.}_{#nu} [GeV]";
    // Add more mappings as needed

    // Fallback: just return var unchanged if no match found
    return var;
}


void Save1DSlices(TH2D* h2d, bool slice_along_q3, const std::string& label, TCanvas* c, const std::string& pdfOut,const std::string& tag,const std::string& xsec_var1,
                  const std::string& xsec_var2) {

  
    int nSlices = slice_along_q3 ? h2d->GetNbinsY() : h2d->GetNbinsX();
    std::string axisTitle = slice_along_q3 ? formatAxisTitle(xsec_var1): formatAxisTitle(xsec_var2);
    

    for (int i = 1; i <= nSlices; ++i) {
        std::string slice_name = Form("%s_slice_%d", h2d->GetName(), i);
        TH1D* h_slice = slice_along_q3 ? h2d->ProjectionX(slice_name.c_str(), i, i)
                                       : h2d->ProjectionY(slice_name.c_str(), i, i);

        // Set axis labels
        h_slice->SetTitle(Form("%s Slice %d", label.c_str(), i));
        h_slice->GetXaxis()->SetTitle(axisTitle.c_str());
        h_slice->GetYaxis()->SetTitle("Posterior Value");

        h_slice->Draw("HIST E");
        c->Print(pdfOut.c_str());
    }
}



int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mcmc_output.root> <bin_config.yaml>" << std::endl;
        return 1;
    }
 
    

    // === MCMC setup ===
  std::vector<BinDef> binDefs;
  std::vector<double> q0_edges;
  std::vector<double> q3_edges;
  
  std::string yaml_file = argv[1];          // YAML config
  std::string mcmc_file = argv[2];          // ROOT file

  manager* FitManager = new manager(yaml_file);  // <-- correct!!

  if (!FitManager) {
    std::cerr << "Error: Failed to create FitManager from YAML config: " << yaml_file << std::endl;
    return 1;
}


  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();
 auto OutputFile = std::unique_ptr<TFile>(TFile::Open( OutputFileName .c_str(), "RECREATE"));

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
 

 

  xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();

  if (xsec_yaml.empty()) {
    std::cerr << "Error: xsec_yaml path is empty or missing in config.\n";
    return 1;
}


  BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);
  binDefs = binning.binDefs;   // <<---- Add this line to populate binDefs

  q0_edges = binning.q0_edges;
  q3_edges = binning.q3_edges;

  if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Error: Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }


  TAxis* axisX = new TAxis(binning.q0_edges.size() - 1, binning.q0_edges.data());
  TAxis* axisY = new TAxis(binning.q3_edges.size() - 1, binning.q3_edges.data());

   
  //////////////////MaCh3 stuff
  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;
  std::vector<samplePDFFDBase*> DUNEPdfs;
 
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);
  std::unique_ptr<mcmc> MaCh3Fitter = std::make_unique<mcmc>(FitManager);
  std::cout << "line 217" << std::endl;

    // Lets draw the event rate distribution in terms of the xsec variables...

    TH2D* h_q0q3_combined = nullptr;              
    TFile* fout = OutputFile.get();     

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
    
    if (h_q0q3_combined) {
        TFile outFile("my_eventratehistogram_templateparams_Enu_Enubias.root", "RECREATE");
        h_q0q3_combined->GetXaxis()->SetTitle(xsec_var1.c_str());
        h_q0q3_combined->GetYaxis()->SetTitle(xsec_var2.c_str());
        h_q0q3_combined->Write("xsec_eventrate_histo");
        outFile.Close();
    } else {
        std::cerr << "Error: h_q0q3_combined is null — no 2D histograms were produced.\n";
        return 1;
    }
  }
   
    
// Load posterior



TFile* f_post = new TFile(mcmc_file.c_str());
TTree* post = (TTree*)f_post->Get("posteriors");
if (!post) {
    std::cerr << "Error: Could not find 'posteriors' TTree in file " << mcmc_file << std::endl;
    return 1;
}

 int nParams = 0;

// Detect number of xsec_* branches
TObjArray* branches = post->GetListOfBranches();
for (int i = 0; i < branches->GetEntries(); ++i) {
    std::string name = branches->At(i)->GetName();
    if (name.find("xsec_") == 0) {
        int idx = std::stoi(name.substr(5));
        if (idx >= nParams)
            nParams = idx + 1;
    }
}
std::cout << "Detected nParams = " << nParams << std::endl;

// Print all branches starting with "xsec_"
//std::cout << "Branches in posteriors TTree starting with 'xsec_':\n";
for (int i = 0; i < branches->GetEntries(); ++i) {
    std::string name = branches->At(i)->GetName();
    if (name.find("xsec_") == 0) {
       
    }
}


// Print all bin indices from binDefs
/*
std::cout << "binDefs indices found:\n";
for (const auto& bin : binDefs) {
    std::cout << "  " << bin.index << "\n";
}*/

std::vector<double> xsec_vals(nParams, 0.0);
std::vector<double*> xsec_ptrs(nParams, nullptr);

for (int i = 0; i < nParams; ++i)
    xsec_ptrs[i] = &xsec_vals[i];

post->SetBranchStatus("*", 0); // Disable all branches first

for (const auto& bin : binDefs) {
    std::string name = "xsec_" + std::to_string(bin.index);
    post->SetBranchStatus(name.c_str(), 1);
    post->SetBranchAddress(name.c_str(), xsec_ptrs[bin.index]); // critical line
}


// Allocate vectors
//std::vector<double> xsec_vals(nParams, 0.0);
//std::vector<double*> xsec_ptrs(nParams, nullptr);
//for (int i = 0; i < nParams; ++i)
  //  xsec_ptrs[i] = &xsec_vals[i];

// Disable all branches
//post->SetBranchStatus("*", 1);

// Optionally: Fill histograms of mean/stddev
TH2D* h_mean = new TH2D("h_param_mean", "Posterior Mean;q_{0} [GeV];q_{3} [GeV]",
                        q0_edges.size()-1, &q0_edges[0],  // X axis: q0
                        q3_edges.size()-1, &q3_edges[0]); // Y axis: q3

TH2D* h_stddev = new TH2D("h_param_stddev", "Posterior StdDev;",
                          q0_edges.size()-1, &q0_edges[0],
                          q3_edges.size()-1, &q3_edges[0]);

std::vector<double> sum(nParams, 0.0);
std::vector<double> sq_sum(nParams, 0.0);

for (const auto& bin : binDefs) {
    if (bin.index >= nParams) {
        std::cerr << "ERROR: bin.index = " << bin.index << " exceeds nParams = " << nParams << std::endl;
    }
}
/*
std::cout << "First 5 entries:\n";
for (int i = 0; i < 5; ++i) {
    post->GetEntry(i);
    std::cout << "Entry " << i << ": ";
    for (int j = 0; j < 5; ++j) {
        std::cout << "xsec_" << j << " = " << xsec_vals[j] << " ";
    }
    std::cout << std::endl;
}*/


Long64_t nEntries = post->GetEntries();
for (Long64_t i = 0; i < nEntries; ++i) {
    post->GetEntry(i);
    for (const auto& bin : binDefs) {
        double val = xsec_vals[bin.index];
        //std::cout << "Entry " << i << ", Bin " << bin.index << ", xsec = " << val << std::endl;
        sum[bin.index] += val;
        sq_sum[bin.index] += val * val;
    }
}
 

for (const auto& bin : binDefs) {
    double nEntries = post->GetEntries() ;
    std::cout << "Bin index: " << bin.index
              << " q0: [" << bin.q0_min << ", " << bin.q0_max << "]"
              << " q3: [" << bin.q3_min << ", " << bin.q3_max << "]\n";
    double mean = sum[bin.index] / nEntries;
    double var = std::max(0.0, (sq_sum[bin.index] / nEntries) - mean * mean);
    double stddev = std::sqrt(var);

    std::cout << "mean = " << mean << "in bin = " <<  bin.index << std::endl;
    double q0_center = 0.5 * (bin.q0_min + bin.q0_max);
    double q3_center = 0.5 * (bin.q3_min + bin.q3_max);

    int bin_q0 = h_mean->GetXaxis()->FindBin(q0_center); // q0 = x
    int bin_q3 = h_mean->GetYaxis()->FindBin(q3_center); // q3 = y

    h_mean->SetBinContent(bin_q0, bin_q3, mean);
    h_stddev->SetBinContent(bin_q0, bin_q3, stddev);
}

// Save histograms
h_mean->Write("xsec_param_mean");
h_stddev->Write("xsec_param_stddev");


   
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetRightMargin(0.15);

    std::string pdfOut = AddTimestampToFilename("mcmc_summary_plots.pdf");
    c->Print((pdfOut + "[").c_str());
    c->Clear(); h_q0q3_combined->Draw("COLZ"); c->Print(pdfOut.c_str());
    c->Clear(); h_mean->Draw("COLZ"); c->Print(pdfOut.c_str());
    c->Clear(); h_stddev->Draw("COLZ"); c->Print(pdfOut.c_str());
   

    // Slice vertically: q3 fixed, q0 varying (i.e., slices along q3 axis)
    Save1DSlices(h_mean, true, "Posterior Mean", c, pdfOut," slices along q3 axis",
             xsec_var1, xsec_var2);
    Save1DSlices(h_stddev, true, "Posterior StdDev", c, pdfOut, "slices along q3 axis",
             xsec_var1, xsec_var2);

    //Also slice horizontally (q0 fixed, q3 varying)
    Save1DSlices(h_mean, false, "Posterior Mean (q0 slices)", c, pdfOut,"slices along q0 axis",
             xsec_var1, xsec_var2);
    Save1DSlices(h_stddev, false, "Posterior StdDev (q0 slices)", c, pdfOut,"slices along q0 axis",
             xsec_var1, xsec_var2);


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
