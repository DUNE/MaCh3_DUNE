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
// ADD THIS
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
if (cut["TrueNeutrinoEnergy"]) {
const auto& q0 = cut["TrueNeutrinoEnergy"];
if (q0.IsSequence() && q0.size() == 2) {
q0_min = q0[0].as<double>();
q0_max = q0[1].as<double>();
q0_set.insert(q0_min);
q0_set.insert(q0_max);
}
}
if (cut["ERecProxy_minus_Enu"]) {
const auto& q3 = cut["ERecProxy_minus_Enu"];
if (q3.IsSequence() && q3.size() == 2) {
q3_min = q3[0].as<double>();
q3_max = q3[1].as<double>();
q3_set.insert(q3_min);
q3_set.insert(q3_max);
}
}
}

/*
if (q0_min >= 0 && q3_min >= 0) {
binDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
}*/


if (std::isfinite(q0_min) && std::isfinite(q3_min)) {
    binDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
}
}

q0_edges.assign(q0_set.begin(), q0_set.end());
q3_edges.assign(q3_set.begin(), q3_set.end());

std::cout << "Parsed " << binDefs.size() << " bins from YAML\n";
std::cout << "   q0 bins: " << q0_edges.size() - 1
<< ", q3 bins: " << q3_edges.size() - 1 << "\n";
} catch (const std::exception& e) {
std::cerr << " YAML parsing failed: " << e.what() << "\n";
}
}



void FixLowStatParams(covarianceXsec* xsec, TH2D* h_q0q3, double total_events, double frac_threshold, std::vector<std::string>& fixed_names_out, const std::vector<BinDef>& binDefs) {
    double threshold = frac_threshold * total_events;
    std::cout << "\n🔍 Threshold for freezing: " << threshold << " events\n";

    for (const auto& bin : binDefs) {
        double q0 = 0.5 * (bin.q0_min + bin.q0_max);
        double q3 = 0.5 * (bin.q3_min + bin.q3_max);
        int bin_q0 = h_q0q3->GetYaxis()->FindBin(q0);
        int bin_q3 = h_q0q3->GetXaxis()->FindBin(q3);

        if (bin_q0 <= 0 || bin_q3 <= 0 || bin_q0 > h_q0q3->GetNbinsY() || bin_q3 > h_q0q3->GetNbinsX()) {
            std::cerr << " Bin index out of bounds: q0_bin=" << bin_q0 << ", q3_bin=" << bin_q3 << std::endl;
            continue;
        }

        double val = h_q0q3->GetBinContent(bin_q3, bin_q0);

        std::cout << "Param " << bin.index << " at (q0=" << q0 << ", q3=" << q3 << ") → Bin content: " << val;
        if (val < threshold) {
            std::cout << " → Frozen";
            xsec->setSingleParameter(bin.index, 0.0);
            if (!xsec->isParameterFixed(bin.index)) xsec->toggleFixParameter(bin.index);
            fixed_names_out.push_back(xsec->GetParFancyName(bin.index));
        }
        std::cout << std::endl;
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

  std::string xsec_param_yaml = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_5.0GeV_fullgrid_smallerbins.yaml";

  std::ifstream test("/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_5.0GeV_fullgrid_smallerbins.yaml");
    if (!test.is_open()) {
    std::cerr << "File not found: /scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_5.0GeV_fullgrid_smallerbins.yaml" << std::endl;
    return 1;
  }
  extract_q0q3_bins_from_yaml(xsec_param_yaml, binDefs, q0_edges, q3_edges);



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
  /*
  for (auto& pdf : DUNEPdfs) {
      //TH2D* h = dynamic_cast<samplePDFDUNEBeamFD*>(pdf)->GetQ0Q3FillingHist();
       auto* typed_pdf = dynamic_cast<samplePDFDUNEBeamFD*>(pdf);
      TH2D* h = typed_pdf->GetQ0Q3FillingHist();
      std::cout << "Hist from PDF integral: " << h->Integral() << std::endl;

      if (!h_q0q3_combined)
          h_q0q3_combined = (TH2D*)h->Clone("h_q0q3_combined");
      else
          h_q0q3_combined->Add(h);
  }
  */

  for (auto& pdf : DUNEPdfs) {
    auto* dunePdf = dynamic_cast<samplePDFDUNEBeamFD*>(pdf);
    TH2D* h = dunePdf->GetQ0Q3FillingHist();
    if (!h_q0q3_combined)
        h_q0q3_combined = (TH2D*)h->Clone("h_q0q3_combined");
    else
        h_q0q3_combined->Add(h);
}

  double h_q0q3_combined_int = h_q0q3_combined->Integral();
  std::cout << "Integral of h_q0q3_combined = " << h_q0q3_combined_int << std::endl;


  // === Step 1: Build flattened prediction ===
  TH1D* h_flat = (TH1D*)DUNEPdfs[0]->get1DHist()->Clone("h_flat");
  for (size_t i = 1; i < DUNEPdfs.size(); ++i)
    h_flat->Add(DUNEPdfs[i]->get1DHist());

  
  std::vector<std::string> fixed_names;

  if (!StartFromPreviousChain) xsec->throwParameters();
  //FixLowStatParams(xsec, h_q0q3, total, 0.0001, fixed_names, binDefs);
  FixLowStatParams(xsec, h_q0q3_combined, h_q0q3_combined->Integral(), 0.0001, fixed_names, binDefs);

  std::cout << "Fixed " << fixed_names.size() << " low-stat parameters." << std::endl;

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
  TH2D* h_event_rate = (TH2D*)h_q0q3_combined->Clone("h_event_rate");
  TH2D* h_param_init  = (TH2D*)h_q0q3_combined->Clone("h_param_init");
  TH2D* h_param_frozen = (TH2D*)h_q0q3_combined->Clone("h_param_frozen");
  TH2D* h_param_mean   = (TH2D*)h_q0q3_combined->Clone("h_param_mean");
  TH2D* h_param_stddev = (TH2D*)h_q0q3_combined->Clone("h_param_stddev");

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
  if (q0_edges.empty() || q3_edges.empty()) {
    std::cerr << " Error: q0 or q3 edges are empty. Histogram cannot be created." << std::endl;
    return 1;
  }
  //std::cout << "[DEBUG] About to clone h_q0q3: " << h_q0q3 << ", title: " << h_q0q3->GetTitle() << std::endl;


  // === Fill init/frozen state ===
  for (const auto& bin : binDefs) {
    double q0 = 0.5 * (bin.q0_min + bin.q0_max);
    double q3 = 0.5 * (bin.q3_min + bin.q3_max);
    int bin_q0 = h_q0q3_combined ->GetYaxis()->FindBin(q0);
    int bin_q3 = h_q0q3_combined ->GetXaxis()->FindBin(q3);

    if (bin_q0 <= 0 || bin_q0 > h_q0q3_combined ->GetNbinsY()) continue;
    if (bin_q3 <= 0 || bin_q3 > h_q0q3_combined ->GetNbinsX()) continue;

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
    int bin_q0 = h_q0q3_combined ->GetYaxis()->FindBin(q0);
    int bin_q3 = h_q0q3_combined ->GetXaxis()->FindBin(q3);

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
  std::string pdfOutName = AddTimestampToFilename("param_summary_q0q3_someparamsfrozen.pdf");
  c->Print((pdfOutName + "[").c_str());  // Open multi-page PDF

  c->Clear(); h_event_rate->Draw("COLZ");     c->Print(pdfOutName.c_str());
  c->Clear(); h_param_init->Draw("COLZ");     c->Print(pdfOutName.c_str());
  c->Clear(); h_param_frozen->Draw("COLZ");   c->Print(pdfOutName.c_str());
  c->Clear(); h_param_mean->Draw("COLZ");     c->Print(pdfOutName.c_str());
  c->Clear(); h_param_stddev->Draw("COLZ");   c->Print(pdfOutName.c_str());

  c->Print((pdfOutName + "]").c_str());  // Close multi-page PDF


  //TFile* fOut = new TFile("param_summary_q0q3_someparamsfrozen.root", "RECREATE");
  std::string rootOutName = AddTimestampToFilename("param_summary_q0q3_someparamsfrozen.root");
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


}
