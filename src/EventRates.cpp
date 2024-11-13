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

#include "samplePDFDUNE/MaCh3DUNEFactory.h"

void WriteHistogramsToFile(std::string OutFileName, std::vector<TH1D*> Histograms, std::vector<TH2D*> Histograms2D) {

  // Now write out the saved hsitograms to file 
  auto OutputFile = std::unique_ptr<TFile>(new TFile(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for(auto Hist : Histograms){
    Hist->Write();
  }
  for(auto Hist : Histograms2D){
    Hist->Write();
  }
  OutputFile->Close();

  return;
}


void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1D*> Histograms) {

  // Now write out the saved hsitograms to file 


  //Remove root from end of file
  OutFileName.erase(OutFileName.find('.'));
  OutFileName+=".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName+"[").c_str());
  for(auto Hist : Histograms){
    Hist->Draw("HIST");
	c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName+"]").c_str());

  return;
}

int main(int argc, char * argv[]) {
  if(argc == 1){
    MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
    return 1;
  }

  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));

  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFFD objects
  
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc); 

  std::vector<TH1D*> DUNEHists;
  std::vector<TH2D*> DUNE2DHists;
  for(auto Sample : DUNEPdfs){
    Sample->reweight();
    if (Sample->GetNDim() == 2){
      DUNE2DHists.push_back(Sample->get2DHist());
      }
    DUNEHists.push_back(Sample->get1DHist());
    std::cout << "Sample LLH: " << Sample->GetName() << " has likelihood " << Sample->GetLikelihood() << std::endl;
    
    std::string EventRateString = fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetName(), EventRateString);
  }

  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRatesOutput.root");

  WriteHistogramsToFile(OutFileName, DUNEHists, DUNE2DHists); 
  Write1DHistogramsToPdf(OutFileName, DUNEHists);
}
