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

struct BinDef {
  int index;
  double q0_min, q0_max;
  double q3_min, q3_max;
};


void extract_xsec_bins_and_edges(
  const std::string& filename,
  std::vector<BinDef>& binDefs, 
  std::vector<double>& q0Edges,
  std::vector<double>& q3Edges)
{
  std::ifstream file(filename);
  std::set<double> q0Set, q3Set;
  binDefs.clear();

  if (!file.is_open()) {
    std::cerr << "Failed to open: " << filename << std::endl;
    return;
  }

  std::string line;
  bool isFirst = true;
  while (std::getline(file, line)) {
    if (isFirst) { isFirst = false; continue; }

    std::stringstream ss(line);
    std::string token;
    BinDef bin;
    int col = 0;

    // Skip the first token (ParameterName)
    std::getline(ss, token, ',');

    while (std::getline(ss, token, ',')) {
      token.erase(0, token.find_first_not_of(" \t"));
      token.erase(token.find_last_not_of(" \t") + 1);
      switch (col++) {
        case 0: bin.q0_min = std::stod(token); break;
        case 1: bin.q0_max = std::stod(token); break;
        case 2: bin.q3_min = std::stod(token); break;
        case 3: bin.q3_max = std::stod(token); break;
      }
      std::cout <<  "In extract_xsec_binedges_fromedges "<<"[DEBUG] q0: " << bin.q0_min << " → " << bin.q0_max<< ", q3: " << bin.q3_min << " → " << bin.q3_max << std::endl;
    }

    binDefs.push_back(bin);
    q0Set.insert(bin.q0_min); q0Set.insert(bin.q0_max);
    q3Set.insert(bin.q3_min); q3Set.insert(bin.q3_max);
  }

  q0Edges.assign(q0Set.begin(), q0Set.end());
  q3Edges.assign(q3Set.begin(), q3Set.end());

}


/*
/// --- Utility: Read q0/q3 edges from txt binning file ---
void extract_xsec_binedges_fromtxt(const std::string& filename, std::vector<double>& q0Edges, std::vector<double>& q3Edges) {
  std::ifstream file(filename);
  std::set<double> q0Set, q3Set;
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  std::string line;
  bool firstLine = true;
  while (std::getline(file, line)) {
    if (firstLine) { firstLine = false; continue; }
    std::stringstream ss(line);
    std::string token;
    int colIndex = 0;
    double q0_min = 0, q0_max = 0, q3_min = 0, q3_max = 0;

    while (std::getline(ss, token, ',')) {
      token.erase(0, token.find_first_not_of(" \t"));
      token.erase(token.find_last_not_of(" \t") + 1);
      switch (colIndex++) {
        case 1: q0_min = std::stod(token); break;
        case 2: q0_max = std::stod(token); break;
        case 3: q3_min = std::stod(token); break;
        case 4: q3_max = std::stod(token); break;
      }
    }
    q0Set.insert(q0_min); q0Set.insert(q0_max);
    q3Set.insert(q3_min); q3Set.insert(q3_max);
  }

  q0Edges.assign(q0Set.begin(), q0Set.end());
  q3Edges.assign(q3Set.begin(), q3Set.end());
}*/

void extract_xsec_binedges_fromtxt(const std::string& filename,std::vector<double>& q0Edges,std::vector<double>& q3Edges) {
  std::ifstream file(filename);
  std::set<double> q0Set, q3Set;

  if (!file.is_open()) {
  std::cerr << "Failed to open file: " << filename << std::endl;
  return;
  }

  std::string line;
  bool firstLine = true;
  while (std::getline(file, line)) {
    if (firstLine) {
    firstLine = false;
    continue;
    }

    std::stringstream ss(line);
    std::string token;
    int colIndex = 0;
    double q0_min = 0, q0_max = 0, q3_min = 0, q3_max = 0;

    try {
      while (std::getline(ss, token, ',')) {
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);

        switch (colIndex++) {
        case 1: q0_min = std::stod(token); break;
        case 2: q0_max = std::stod(token); break;
        case 3: q3_min = std::stod(token); break;
        case 4: q3_max = std::stod(token); break;
        }
      }
    } catch (const std::invalid_argument& e) {
  std::cerr << "Invalid number in line: " << line << std::endl;
  continue;  // Skip malformed line
  }

  q0Set.insert(q0_min);
  q0Set.insert(q0_max);
  q3Set.insert(q3_min);
  q3Set.insert(q3_max);
  std::cout <<  "In extract_xsec_binedges_fromtxt "<<"[DEBUG] q0: " << q0_min << " → " << q0_max << ", q3: " << q3_min << " → " << q3_max << std::endl;
  }

  q0Edges.assign(q0Set.begin(), q0Set.end());
  q3Edges.assign(q3Set.begin(), q3Set.end());
  

}


/// --- Utility: Freeze low-stat parameters based on bin content ---
void FixLowStatParams(covarianceXsec* xsec, TH2D* h_q0q3, double total_events, double frac_threshold, std::vector<std::string>& fixed_names_out) {
  int nbins_q0 = h_q0q3->GetNbinsY();
  int nbins_q3 = h_q0q3->GetNbinsX();
  int nParams  = xsec->GetNumParams();
  double threshold = frac_threshold * total_events;

  for (int i = 0; i < nParams; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;
    if (bin_q0 > nbins_q0 || bin_q3 > nbins_q3) continue;

    double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
    if (q0 > q3) continue;

    double val = h_q0q3->GetBinContent(bin_q3, bin_q0);
    if (val < threshold) {
      xsec->setSingleParameter(i, 0.0);
      xsec->toggleFixParameter(i);
      fixed_names_out.push_back(xsec->GetParFancyName(i));
    }
  }
}

/// --- Main Fit and Diagnostics ---
int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: Fit config.yaml" << std::endl;
    return 1;
  }

  // === MCMC setup ===
  
  std::vector<BinDef> binDefs;
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
  std::vector<double> q0_edges, q3_edges;
  std::string xsec_binning_txtfile = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/utils/q0q3_0.0_5.0GeV_triangle_parameter_list_morebins.txt";
  extract_xsec_binedges_fromtxt(xsec_binning_txtfile, q3_edges, q0_edges);

  std::cout << "[DEBUG] q0_edges.size() = " << q0_edges.size() << std::endl;
  std::cout << "[DEBUG] q3_edges.size() = " << q3_edges.size() << std::endl;


  TH2D* h_q0q3 = new TH2D("h_q0q3", "q_{0} vs q_{3}; q_{3} [GeV]; q_{0} [GeV]",
                          q3_edges.size() - 1, &q3_edges[0],
                          q0_edges.size() - 1, &q0_edges[0]);

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
  FixLowStatParams(xsec, h_q0q3, total, 0.0001, fixed_names);
  
  extract_xsec_bins_and_edges("/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/utils/q0q3_0.0_5.0GeV_triangle_parameter_list_morebins.txt", binDefs, q0_edges, q3_edges);

  // === Histograms ===
  //int nParams = xsec->GetNumParams();
  
  //std::cout << "[DEBUG] nbins: " << h_param_init->GetNbinsX() << "x" << h_param_init->GetNbinsY() << std::endl;



  // === Record frozen parameters ===
  std::ofstream fixedOut("fixed_parameters.txt");
  for (const auto& name : fixed_names) fixedOut << name << "\n";
  fixedOut.close();

  // === Record initial and frozen parameter values ===
  /*for (const auto& bin : binDefs) {
    double q0 = 0.5 * (bin.q0_min + bin.q0_max);
    double q3 = 0.5 * (bin.q3_min + bin.q3_max);
    int bin_q0 = h_q0q3->GetYaxis()->FindBin(q0);
    int bin_q3 = h_q0q3->GetXaxis()->FindBin(q3);

    if (bin_q0 <= 0 || bin_q0 > h_q0q3->GetNbinsY()) continue;
    if (bin_q3 <= 0 || bin_q3 > h_q0q3->GetNbinsX()) continue;

    // All parameters start at 1.0 (before freezing)
    //h_param_init->SetBinContent(bin_q3, bin_q0, 1.0);

    // Frozen ones get set to 0.0
    //double frozen_val = xsec->isParameterFixed(bin.index) ? 0.0 : 1.0;
    //h_param_frozen->SetBinContent(bin_q3, bin_q0, frozen_val);
  }*/
  //std::cout << "[CHECK] h_param_frozen integral: " << h_param_frozen->Integral() << std::endl;
  //std::cout << "[CHECK] h_param_init integral: " << h_param_init->Integral() << std::endl;


  //h_param_frozen->SetMinimum(0.0);
  //h_param_frozen->SetMaximum(1.0);

  //h_param_init->SetMinimum(0.0);
  //h_param_init->SetMaximum(1.0);

  
  
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

  MaCh3Fitter->runMCMC();


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


// === Prepare Histograms ===
TH2D* h_event_rate = (TH2D*)h_q0q3->Clone("h_event_rate");
TH2D* h_param_init  = (TH2D*)h_q0q3->Clone("h_param_init");
TH2D* h_param_frozen = (TH2D*)h_q0q3->Clone("h_param_frozen");
TH2D* h_param_mean   = (TH2D*)h_q0q3->Clone("h_param_mean");
TH2D* h_param_stddev = (TH2D*)h_q0q3->Clone("h_param_stddev");

h_param_init->Reset();
h_param_frozen->Reset();
h_param_mean->Reset();
h_param_stddev->Reset();

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

std::cout << "\n✅ Saved:\n - param_summary_q0q3.pdf\n - param_summary_q0q3.root\n" << std::endl;

}
