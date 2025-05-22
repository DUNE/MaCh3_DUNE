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

#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>

void extract_xsec_binedges_fromtxt(const std::string& filename, std::vector<double>& q0Edges,std::vector<double>& q3Edges) {
  std::ifstream file(filename);
  std::set<double> q0Set;
  std::set<double> q3Set;

  if (!file.is_open()) {
  std::cerr << "Failed to open file: " << filename << std::endl;
  return;
  }

  std::string line;
  bool isFirstLine = true;
  while (std::getline(file, line)) {
  if (isFirstLine) {
  isFirstLine = false;
  continue; // Skip header
  }

  std::stringstream ss(line);
  std::string token;
  int colIndex = 0;
  double q0_min = 0.0, q0_max = 0.0;
  double q3_min = 0.0, q3_max = 0.0;

  while (std::getline(ss, token, ',')) {
  // Trim whitespace
  token.erase(0, token.find_first_not_of(" \t"));
  token.erase(token.find_last_not_of(" \t") + 1);

  switch (colIndex) {
  case 1: q0_min = std::stod(token); break;
  case 2: q0_max = std::stod(token); break;
  case 3: q3_min = std::stod(token); break;
  case 4: q3_max = std::stod(token); break;
  }
  ++colIndex;
  }

  q0Set.insert(q0_min);
  q0Set.insert(q0_max);
  q3Set.insert(q3_min);
  q3Set.insert(q3_max);
  }

  q0Edges.assign(q0Set.begin(), q0Set.end());
  q3Edges.assign(q3Set.begin(), q3Set.end());
}

void FixLowStatParams(covarianceXsec* xsec,
  TH2D* h_q0q3,
  double total_events,
  double frac_threshold,
  std::vector<std::string>& fixed_names_out) {
int nbins_q0 = h_q0q3->GetNbinsY();
int nbins_q3 = h_q0q3->GetNbinsX();
int nParams  = xsec->GetNumParams();
double threshold = frac_threshold * total_events;

int nFixed = 0;

for (int i = 0; i < nParams; ++i) {
int bin_q0 = i / nbins_q3 + 1;
int bin_q3 = i % nbins_q3 + 1;

if (bin_q0 > nbins_q0 || bin_q3 > nbins_q3) continue;

double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);

if (q0 > q3) continue;  // unphysical

double bin_val = h_q0q3->GetBinContent(bin_q3, bin_q0);
if (bin_val < threshold) {
xsec->setSingleParameter(i, 0.0);
xsec->toggleFixParameter(i);

++nFixed;
std::string name = xsec->GetParFancyName(i);
fixed_names_out.push_back(name);

MACH3LOG_INFO("FIXED param '{}' (index {}): bin content {:.4f} < {:.4f}",
name, i, bin_val, threshold);
}
}

MACH3LOG_INFO("Fixed {} of {} xsec parameters (bins < {:.2f}% of total)",
nFixed, nParams, frac_threshold * 100.0);
}




int main(int argc, char * argv[]) {

  auto FitManager = MaCh3ManagerFactory(argc, argv);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  ParameterHandlerGeneric* xsec = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  //Some place to store the histograms
  std::vector<TH1*> PredictionHistograms;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {

    std::string name = DUNEPdfs[sample_i]->GetTitle();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());

    DUNEPdfs[sample_i] -> Reweight();
    PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i]->GetMCHist(DUNEPdfs[sample_i]->GetNDim())->Clone(NameTString+"_unosc")));

    if (DUNEPdfs[sample_i]->GetNDim() == 1){
      DUNEPdfs[sample_i]->AddData(static_cast<TH1D*>(PredictionHistograms[sample_i]));
    } else if (DUNEPdfs[sample_i]->GetNDim() == 2){
      DUNEPdfs[sample_i]->AddData(static_cast<TH2D*>(PredictionHistograms[sample_i]));
    }

    else {
      MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting");
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }


  }
  // Step 0 -1D flatnenned histo
  TH1D* h_flat = (TH1D*)DUNEPdfs[0]->get1DHist()->Clone("h_flat");
  for (size_t i = 1; i < DUNEPdfs.size(); ++i) {
    h_flat->Add(DUNEPdfs[i]->get1DHist());
  }
  ////////////////// First of all, we need the xsec bin edges - where are each of the bins in the xsec model?
  std::string xsec_binning_txtfile = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/utils/q0q3_0.0_5.0GeV_triangle_parameter_list_morebins.txt";
  std::vector<double> q0_edges;
  std::vector<double> q3_edges;
  extract_xsec_binedges_fromtxt(xsec_binning_txtfile, q3_edges, q0_edges);

   // Now make a TH2D that matches the binning
   TH2D* h_q0q3 = new TH2D("h_q0q3", "q_{0} vs q_{3};q_{3} [GeV];q_{0} [GeV]",
    q3_edges.size() - 1, &q3_edges[0],
    q0_edges.size() - 1, &q0_edges[0]);

  /// Now map
  int nbins_q3 = q3_edges.size() - 1;
  int nbins_q0 = q0_edges.size() - 1;

  for (int i = 0; i < nbins_q0 * nbins_q3; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;

    double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);

    if (q0 <= q3) {
      h_q0q3->SetBinContent(bin_q3, bin_q0, h_flat->GetBinContent(i + 1));
    } else {
      h_q0q3->SetBinContent(bin_q3, bin_q0, 0.0);  // mask unphysical bins
    }
  }
  double total = h_q0q3->Integral();
  std::cout << "total, h_q0q3->Integral() = " << total << std::endl;

  std::vector<std::string> fixed_names;
  FixLowStatParams(xsec, h_q0q3, total, 0.001, fixed_names);

  std::ofstream fixedOut("fixed_parameters.txt");
  for (const auto& name : fixed_names) fixedOut << name << "\n";
  fixedOut.close();




/*
// Optional: replace bin labels with parameter names (if they match)
for (int i = 0; i < xsec->GetNumParams(); ++i) {
  hNominal->GetXaxis()->SetBinLabel(i+1, xsec->GetParFancyName(i).c_str());
}
hNominal->SetTitle("Nominal Prediction per Parameter Bin;Parameter;Events");
hNominal->SetStats(0);
TCanvas* c = new TCanvas("c", "Nominal Event Rate", 1200, 600);
hNominal->Draw("HIST");
c->SaveAs("param_bin_events.pdf");*/


  //###########################################################################################################
  //MCMC

  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);
  //Start chain from random position unless continuing a chain
  if(!StartFromPreviousChain){
    if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false)) {
      xsec->ThrowParameters();
    }
  }

    // Now reset fixed parameters to nominal AFTER throw
  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    if (std::find(fixed_names.begin(), fixed_names.end(), xsec->GetParFancyName(i)) != fixed_names.end()) {
      xsec->setSingleParameter(i, 0.0);
      if (!xsec->isParameterFixed(i)) xsec->toggleFixParameter(i);
      MACH3LOG_INFO("Re-fixed '{}' after throw", xsec->GetParFancyName(i));
    }
  }

  //Add systematic objects
  MaCh3Fitter->AddSystObj(xsec);

  if (StartFromPreviousChain) {
    std::string PreviousChainPath = FitManager->raw()["General"]["PosFileName"].as<std::string>();
    MACH3LOG_INFO("MCMC getting starting position from: {}",PreviousChainPath);
    MaCh3Fitter->StartFromPreviousFit(PreviousChainPath);
  }

  //Add samples
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }
  /*
  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    if (xsec->isParameterFixed(i)) {
      std::cout << "Param " << i << " = " << xsec->CurrVal(i)
                << " (fixed)" << std::endl;
    }
  }

   std::ofstream fout("q0q3_bin_counts.txt");
   fout << "# q0 [GeV], q3 [GeV], Events\n";

   for (int iy = 1; iy <= h_q0q3->GetNbinsY(); ++iy) {
     double q0 = h_q0q3->GetYaxis()->GetBinCenter(iy);
     for (int ix = 1; ix <= h_q0q3->GetNbinsX(); ++ix) {
       double q3 = h_q0q3->GetXaxis()->GetBinCenter(ix);
       double val = h_q0q3->GetBinContent(ix, iy);
       fout << std::fixed << std::setprecision(4) << q0 << ", " << q3 << ", " << val << "\n";
     }
   }
   fout.close();
   MACH3LOG_INFO("Wrote q0-q3 bin event counts to q0q3_bin_counts.txt");
   */

  //Run fit
  MaCh3Fitter->RunMCMC();

  ////////////////////////Make some plots
  // Get binning
  //int nbins_q0 = h_q0q3->GetNbinsY();
  //int nbins_q3 = h_q0q3->GetNbinsX();
  int nParams = xsec->GetNumParams();

  // === Histogram: 1. Event Rate (already created as h_q0q3) ===
  TH2D* h_event_rate = (TH2D*)h_q0q3->Clone("h_event_rate");
  h_event_rate->SetTitle("1. Nominal Event Rate; q_{3} [GeV]; q_{0} [GeV]");

  // === Histogram: 2. Initial Parameter Values (before freezing) ===
  TH2D* h_param_init = (TH2D*)h_q0q3->Clone("h_param_init");
  h_param_init->Reset();
  h_param_init->SetTitle("2. Parameter Values Before Freezing");

  for (int i = 0; i < nParams; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;
    if (bin_q0 > nbins_q0) continue;

    double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
    if (q0 > q3) continue;

    h_param_init->SetBinContent(bin_q3, bin_q0, xsec->getParCurr(i));
  }

  // === Histogram: 3. Parameter Values After Freezing ===
  TH2D* h_param_frozen = (TH2D*)h_q0q3->Clone("h_param_frozen");
  h_param_frozen->Reset();
  h_param_frozen->SetTitle("3. Parameter Values After Freezing");

  for (int i = 0; i < nParams; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;
    if (bin_q0 > nbins_q0) continue;

    double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
    if (q0 > q3) continue;

    h_param_frozen->SetBinContent(bin_q3, bin_q0, xsec->getParCurr(i));
  }

  TFile* chainFile = new TFile(OutputFileName.c_str(), "READ");
  TTree* post = (TTree*)chainFile->Get("posteriors");
  if (!post) {
    std::cerr << "❌ ERROR: Could not find 'posteriors' TTree in " << OutputFileName << std::endl;
    return 1;
  }


  // === Histogram: 4. Posterior Mean ===

  TH2D* h_param_mean = (TH2D*)h_q0q3->Clone("h_param_mean");
  h_param_mean->Reset();
  h_param_mean->SetTitle("4. Posterior Mean Parameter Values");

  std::vector<double> sum(nParams, 0.0);
  int nEntries = post->GetEntries();
  std::vector<double> xsec_vals(nParams, 0.0);
  for (int i = 0; i < nParams; ++i)
    post->SetBranchAddress(("xsec_" + std::to_string(i)).c_str(), &xsec_vals[i]);

  for (int entry = 0; entry < nEntries; ++entry) {
    post->GetEntry(entry);
    for (int i = 0; i < nParams; ++i) sum[i] += xsec_vals[i];
  }

  for (int i = 0; i < nParams; ++i) {
    int bin_q0 = i / nbins_q3 + 1;
    int bin_q3 = i % nbins_q3 + 1;
    if (bin_q0 > nbins_q0) continue;

    double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
    double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
    if (q0 > q3) continue;

    double mean_val = sum[i] / nEntries;
    h_param_mean->SetBinContent(bin_q3, bin_q0, mean_val);
  }

  TCanvas* c = new TCanvas("c_summary", "", 900, 700);

  c->Print("param_summary_q0q3.pdf[");  // open PDF

  h_event_rate->Draw("COLZ");
  c->Print("param_summary_q0q3.pdf");

  h_param_init->Draw("COLZ");
  c->Print("param_summary_q0q3.pdf");

  h_param_frozen->Draw("COLZ");
  c->Print("param_summary_q0q3.pdf");

  h_param_mean->Draw("COLZ");
  c->Print("param_summary_q0q3.pdf");

  c->Print("param_summary_q0q3.pdf]");  // close PDF

  TFile* fout = new TFile("param_summary_q0q3.root", "RECREATE");
  h_event_rate->Write();
  h_param_init->Write();
  h_param_frozen->Write();
  h_param_mean->Write();
  fout->Close();



  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
