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

#include <string>
#include <yaml-cpp/yaml.h>

// 👇 ADD THIS
struct BinDef {
  int index;
  double q0_min, q0_max;
  double q3_min, q3_max;
};

void extract_q0q3_bins_from_yaml(const std::string& yaml_file,
  std::vector<BinDef>& binDefs,
  std::vector<double>& q0_edges,
  std::vector<double>& q3_edges) {
std::set<double> q0_set, q3_set;
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
if (cut["q0"]) {
const auto& q0 = cut["q0"];
if (q0.IsSequence() && q0.size() == 2) {
q0_min = q0[0].as<double>();
q0_max = q0[1].as<double>();
q0_set.insert(q0_min);
q0_set.insert(q0_max);
}
}
if (cut["q3"]) {
const auto& q3 = cut["q3"];
if (q3.IsSequence() && q3.size() == 2) {
q3_min = q3[0].as<double>();
q3_max = q3[1].as<double>();
q3_set.insert(q3_min);
q3_set.insert(q3_max);
}
}
}

if (q0_min >= 0 && q3_min >= 0) {
binDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
}
}

q0_edges.assign(q0_set.begin(), q0_set.end());
q3_edges.assign(q3_set.begin(), q3_set.end());

std::cout << "✅ Parsed " << binDefs.size() << " bins from YAML\n";
std::cout << "   q0 bins: " << q0_edges.size() - 1
<< ", q3 bins: " << q3_edges.size() - 1 << "\n";
} catch (const std::exception& e) {
std::cerr << "❌ YAML parsing failed: " << e.what() << "\n";
}
}


/// --- Utility: Freeze low-stat parameters based on bin content ---
void FixLowStatParams(covarianceXsec* xsec, TH2D* h_q0q3, double total_events, double frac_threshold, std::vector<std::string>& fixed_names_out, const std::vector<BinDef>& binDefs) {
  int nbins_q0 = h_q0q3->GetNbinsY();
  int nbins_q3 = h_q0q3->GetNbinsX();
  int nParams  = xsec->GetNumParams();
  double threshold = frac_threshold * total_events;

  /*
  for (int i = 0; i < nParams; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;
    if (bin_q0 > nbins_q0 || bin_q3 > nbins_q3) continue;

    double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
    if (q0 > q3) continue;

    double val = h_q0q3->GetBinContent(bin_q3, bin_q0);
    std::cout << "Content of bin = " << val << "q0 = " << bin_q0 << "q3 = " << bin_q3 << std::endl;
    if (val < threshold) {
      xsec->setSingleParameter(i, 0.0);
      xsec->toggleFixParameter(i);
      fixed_names_out.push_back(xsec->GetParFancyName(i));
    }
  }*/
  for (const auto& bin : binDefs) {
    double q0 = 0.5 * (bin.q0_min + bin.q0_max);
    double q3 = 0.5 * (bin.q3_min + bin.q3_max);

    int bin_q0 = h_q0q3->GetYaxis()->FindBin(q0);
    int bin_q3 = h_q0q3->GetXaxis()->FindBin(q3);
    double val = h_q0q3->GetBinContent(bin_q3, bin_q0);

    if (val < threshold) {
        xsec->setSingleParameter(bin.index, 0.0);
        xsec->toggleFixParameter(bin.index);
        fixed_names_out.push_back(xsec->GetParFancyName(bin.index));
    }
  }

  
}

/// --- Main Fit and Diagnostics ---
int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: Fit config.yaml" << std::endl;
    return 1;
  }
  TH1::AddDirectory(kFALSE);

  // === MCMC setup ===
  std::vector<BinDef> binDefs;
  std::vector<double> q0_edges, q3_edges;
  //std::vector<BinDef> binDefs;
  manager* FitManager = new manager(argv[1]);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  
  if (!OutputFile || OutputFile->IsZombie()) {
    std::cerr << "Failed to open output file: " << OutputFileName << std::endl;
    return 1;
  }
  OutputFile->cd();


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

  // === Step 1: Build flattened prediction ===
  TH1D* h_flat = (TH1D*)DUNEPdfs[0]->get1DHist()->Clone("h_flat");
  for (size_t i = 1; i < DUNEPdfs.size(); ++i)
    h_flat->Add(DUNEPdfs[i]->get1DHist());

  // === Step 2: Create 2D histogram of (q0,q3) with physical binning ===
  
  std::string xsec_param_yaml = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/q0q3_triangle_1275params.yaml";

  std::ifstream test("/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/q0q3_triangle_1275params.yaml");
    if (!test.is_open()) {
    std::cerr << "❌ File not found: /scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/q0q3_triangle_1275params.yaml" << std::endl;
    return 1;
  }
  extract_q0q3_bins_from_yaml(xsec_param_yaml, binDefs, q0_edges, q3_edges);

 
std::cout << "✅ Extracted " << q0_edges.size() << " q0 edges and " << q3_edges.size() << " q3 edges" << std::endl;


  //std::vector<double> q0_edges, q3_edges;
  
  //extract_q0q3_bins_from_yaml(xsec_param_yaml, q0_edges, q3_edges);

  //extract_q0q3_bins_from_yaml(argv[1], binDefs, q0_edges, q3_edges);

  std::cout << "[DEBUG] q0_edges.size() = " << q0_edges.size() << std::endl;
  std::cout << "[DEBUG] q3_edges.size() = " << q3_edges.size() << std::endl;


  TH2D* h_q0q3 = new TH2D("h_q0q3", "q_{0} vs q_{3}; q_{3} [GeV]; q_{0} [GeV]",
                          q3_edges.size() - 1, &q3_edges[0],
                          q0_edges.size() - 1, &q0_edges[0]);
  h_q0q3->SetDirectory(nullptr);

  int nbins_q3 = q3_edges.size() - 1;
  int nbins_q0 = q0_edges.size() - 1;
  for (int i = 0; i < nbins_q0 * nbins_q3; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;
    double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
    if (q0 <= q3) h_q0q3->SetBinContent(bin_q3, bin_q0, h_flat->GetBinContent(i + 1));
  }

  double total = h_q0q3->Integral();
  std::vector<std::string> fixed_names;

  if (!StartFromPreviousChain) xsec->throwParameters();
  FixLowStatParams(xsec, h_q0q3, total, 0.0001, fixed_names, binDefs);
  
  std::cout << "✅ Fixed " << fixed_names.size() << " low-stat parameters." << std::endl;

  // === Record frozen parameters ===
  std::ofstream fixedOut("fixed_parameters.txt");
  for (const auto& name : fixed_names) fixedOut << name << "\n";
  fixedOut.close();


  // Re-freeze after throw
  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    if (std::find(fixed_names.begin(), fixed_names.end(), xsec->GetParFancyName(i)) != fixed_names.end()) {
      xsec->setSingleParameter(i, 0.0);
      if (!xsec->isParameterFixed(i)) xsec->toggleFixParameter(i);
    }
  }

  if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false))
    MaCh3Fitter->addSystObj(xsec);

  for (auto& pdf : DUNEPdfs) MaCh3Fitter->addSamplePDF(pdf);
  if (StartFromPreviousChain)
    MaCh3Fitter->StartFromPreviousFit(FitManager->raw()["General"]["PosFileName"].as<std::string>());

  // === Prepare Histograms ===
  TH2D* h_event_rate = (TH2D*)h_q0q3->Clone("h_event_rate");
  TH2D* h_param_init  = (TH2D*)h_q0q3->Clone("h_param_init");
  TH2D* h_param_frozen = (TH2D*)h_q0q3->Clone("h_param_frozen");
  TH2D* h_param_mean   = (TH2D*)h_q0q3->Clone("h_param_mean");
  TH2D* h_param_stddev = (TH2D*)h_q0q3->Clone("h_param_stddev");

  h_event_rate->SetDirectory(nullptr);
  h_param_init->SetDirectory(nullptr);
  h_param_frozen->SetDirectory(nullptr);
  h_param_mean->SetDirectory(nullptr);
  h_param_stddev->SetDirectory(nullptr);

  std::cout << "[DEBUG] h_q0q3->IsOnHeap(): " << h_q0q3->IsOnHeap() << std::endl;
  std::cout << "[DEBUG] gDirectory: " << gDirectory->GetName() << std::endl;



  MaCh3Fitter->runMCMC();

  std::cout << "[DEBUG] h_q0q3->IsOnHeap(): " << h_q0q3->IsOnHeap() << std::endl;
  std::cout << "[DEBUG] gDirectory: " << gDirectory->GetName() << std::endl;

  std::cout << "done run MCMC" << std::endl;
  // === Load Posterior ===
  TFile* f_chain = new TFile(OutputFileName.c_str());
  TTree* post = (TTree*)f_chain->Get("posteriors");
  if (!post) {
    std::cerr << "❌ Could not find 'posteriors' TTree in " << OutputFileName << std::endl;
    return 1;
  }
  if (!h_q0q3) {
    std::cerr << "❌ h_q0q3 is null before cloning!" << std::endl;
    return 1;
  }
  if (q0_edges.empty() || q3_edges.empty()) {
    std::cerr << "❌ Error: q0 or q3 edges are empty. Histogram cannot be created." << std::endl;
    return 1;
  }
  //std::cout << "[DEBUG] About to clone h_q0q3: " << h_q0q3 << ", title: " << h_q0q3->GetTitle() << std::endl;


  // === Fill init/frozen state ===
  for (const auto& bin : binDefs) {
    double q0 = 0.5 * (bin.q0_min + bin.q0_max);
    double q3 = 0.5 * (bin.q3_min + bin.q3_max);
    int bin_q0 = h_q0q3->GetYaxis()->FindBin(q0);
    int bin_q3 = h_q0q3->GetXaxis()->FindBin(q3);

    if (bin_q0 <= 0 || bin_q0 > h_q0q3->GetNbinsY()) continue;
    if (bin_q3 <= 0 || bin_q3 > h_q0q3->GetNbinsX()) continue;

    h_param_init->SetBinContent(bin_q3, bin_q0, 1.0);
    double frozen_val = xsec->isParameterFixed(bin.index) ? 0.0 : 1.0;
    h_param_frozen->SetBinContent(bin_q3, bin_q0, frozen_val);
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
    double q3 = 0.5 * (bin.q3_min + bin.q3_max);
    int bin_q0 = h_q0q3->GetYaxis()->FindBin(q0);
    int bin_q3 = h_q0q3->GetXaxis()->FindBin(q3);

    double mean = sum[bin.index] / nEntries;
    double var  = std::max(0.0, (sq_sum[bin.index] / nEntries) - mean * mean);
    double stddev = std::sqrt(var);

    h_param_mean->SetBinContent(bin_q3, bin_q0, mean);
    h_param_stddev->SetBinContent(bin_q3, bin_q0, stddev);
  }

  // === Output diagnostics ===
  std::cout << "[CHECK] h_param_frozen integral: " << h_param_frozen->Integral() << std::endl;
  std::cout << "[CHECK] h_param_init integral: " << h_param_init->Integral() << std::endl;

  // === Save plots ===
  gStyle->SetPalette(kRainBow);
  TCanvas* c = new TCanvas("c", "", 900, 700);
  c->Print("param_summary_q0q3.pdf[");
  h_event_rate->Draw("COLZ");     c->Print("param_summary_q0q3.pdf");
  h_param_init->Draw("COLZ");     c->Print("param_summary_q0q3.pdf");
  h_param_frozen->Draw("COLZ");   c->Print("param_summary_q0q3.pdf");
  h_param_mean->Draw("COLZ");     c->Print("param_summary_q0q3.pdf");
  h_param_stddev->Draw("COLZ");   c->Print("param_summary_q0q3.pdf");
  c->Print("param_summary_q0q3.pdf]");

  TFile* fOut = new TFile("param_summary_q0q3.root", "RECREATE");
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

  std::cout << "\n✅ Saved:\n - param_summary_q0q3.pdf\n - param_summary_q0q3.root\n" << std::endl;

}
