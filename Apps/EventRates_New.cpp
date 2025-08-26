#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TColor.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRint.h>
#include <TStyle.h>

#include "samplePDF/GenericBinningTools.h"

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "mcmc/mcmc.h"

void Write1DHistogramsToFile(std::string OutFileName,
                             std::vector<TH1D *> Histograms) {

  // Now write out the saved hsitograms to file
  auto OutputFile =
      std::unique_ptr<TFile>(new TFile(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for (auto Hist : Histograms) {
    Hist->Write();
  }
  OutputFile->Close();

  return;
}

void Write1DHistogramsToPdf(std::string OutFileName,
                            std::vector<TH1D *> Histograms) {

  // Now write out the saved hsitograms to file

  // Remove root from end of file
  OutFileName.erase(OutFileName.find('.'));
  OutFileName += ".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName + "[").c_str());
  for (auto Hist : Histograms) {
    Hist->Draw("HIST");
    c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName + "]").c_str());

  return;
}


int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "Usage: bin/EventRatesDUNEBeam config.cfg" << std::endl;
    return 1;
  }
  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));
  
  std::string OutFileName = GetFromManager<std::string>(
    fitMan->raw()["General"]["OutputFile"], "EventRates.root");

  // Replace ".root" with "_prism.pdf"
  std::string PrismFileName = OutFileName;
  size_t pos = PrismFileName.find(".root");
  if (pos != std::string::npos) {
      PrismFileName.replace(pos, 5, "_prism.pdf");
  } else {
      PrismFileName += "_prism.pdf";
  }


  covarianceXsec *xsec = nullptr;
  covarianceOsc *osc = nullptr;

 // std::string hists_file (OutFileName + "event_histograms").c_str();
  // ####################################################################################
  // Create samplePDFFD objects

  std::vector<samplePDFFDBase *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  auto gc1 = std::unique_ptr<TCanvas>(new TCanvas("gc1", "gc1", 800, 600));
  gStyle->SetOptStat(false);
  gc1->Print((PrismFileName + "[").c_str());  // Open multi-page PDF
  gc1->Print("GenericBinTest.pdf[");
  

  /*
  if (Sample -> GetNDim() == 1) {
      TH1D *Asimov_1D = (TH1D*)Sample->get1DHist()->Clone(NameTString+"_asimov");
      std::cout << name.c_str() << ": " << Asimov_1D->Integral() << std::endl;
      Sample -> addData(Asimov_1D); 
	}*/

  std::vector<TH1D *> DUNEHists;
  for (auto Sample : DUNEPdfs) {

    // inside the loop over samples
    xsec->setParameters();
    Sample->reweight(); 
    TH1D *Asimov_1D = (TH1D*)Sample->get1DHist()->Clone((Sample->GetTitle()+"_asimov").c_str());
    //Sample->addData(Sample->get1DHist()->Clone((Sample->GetTitle() + "_asimovdata").c_str()));
    std::cout << Sample->GetTitle() 
              << " Asimov integral = " << Asimov_1D->Integral() 
              << " entries = " << Asimov_1D->GetEntries() << std::endl;


    Sample -> addData(Asimov_1D); 
    Sample->reweight();
    xsec->setParameters();
    double nominal =xsec->getNominal(0); //get central value of parameter
    double error = xsec->getDiagonalError(0);
    std::cout<<"nominal  = " << nominal << std::endl; 
    std::cout<<"error  = " << error << std::endl; 
    xsec->setParCurrProp(0, nominal);////////// set //+(2*error)
    std::cout<< "nominal value = " << nominal <<std::endl;;
    double current_value = xsec->getParProp(0);
    std::cout<<"current value  = " << current_value << std::endl; 
    Sample->reweight();

      /*
     if (Sample->GetNDim() == 1) {
       auto myhist2 = (TH1D*)Sample->get1DHist()->Clone((Sample->GetTitle()+"_draw").c_str());
        myhist2->Draw("HIST");
        gc1->Print("GenericBinTest.pdf");
        gc1->Print(PrismFileName.c_str());
        DUNEHists.push_back(myhist2);
      }*/


    if (Sample->generic_binning.GetNDimensions()) {

      auto myhist = GetGenericBinningTH1(*Sample, "myhist");
      myhist->Scale(1, "WIDTH");
      myhist->Draw();
      gc1->Print("GenericBinTest.pdf");
      gc1->Print(PrismFileName.c_str());

     
      if (Sample->generic_binning.GetNDimensions() == 2) {
        auto myhist2 = GetGenericBinningTH2(*Sample, "myhist2");
        myhist2->Draw("COLZ");
        gc1->Print("GenericBinTest.pdf");
        gc1->Print(PrismFileName.c_str());

        for (auto &slice :
             GetGenericBinningTH1Slices(*Sample, 0, "myslicehist")) {
          slice->Draw("colz");
          gc1->Print("GenericBinTest.pdf");
          gc1->Print(PrismFileName.c_str());
        }
      }
      if (Sample->generic_binning.GetNDimensions() == 3) {
        for (auto &slice :
             GetGenericBinningTH2Slices(*Sample, {0, 1}, "myslicehist")) {
          slice->Draw("colz ");
          gc1->Print("GenericBinTest.pdf");
          gc1->Print(PrismFileName.c_str());
        }
      }
    }

    DUNEHists.push_back(Sample->get1DHist());

    std::string EventRateString =
        fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetTitle(),
                  EventRateString);

    std::string LLHString =
        fmt::format("{:.10f}", Sample-> GetLikelihood());
    MACH3LOG_INFO("LLH for {} : {:<5}", Sample->GetTitle(),
                  LLHString);
  }

  gc1->Print("GenericBinTest.root]");
  gc1->Print("GenericBinTest.pdf]");
  gc1->Print((PrismFileName + "]").c_str());

  Write1DHistogramsToFile(OutFileName, DUNEHists);
  Write1DHistogramsToPdf(OutFileName, DUNEHists);
}