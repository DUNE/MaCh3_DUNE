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

void WriteHistogramsToFile(const std::string &OutFileName,
                           const std::vector<TH1 *> &Histograms) {

  std::cout << "[INFO] Writing " << Histograms.size()
            << " histograms to ROOT file: " << OutFileName << std::endl;

  TFile OutputFile(OutFileName.c_str(), "RECREATE");
  if (!OutputFile.IsOpen()) {
    std::cerr << "[ERROR] Could not open " << OutFileName << " for writing.\n";
    return;
  }

  for (auto *Hist : Histograms) {
    if (!Hist) continue;
    std::cout << "  -> Writing " << Hist->GetName() << " (" << Hist->ClassName() << ")\n";
    Hist->Write();
  }

  OutputFile.Write();
  OutputFile.Close();
  std::cout << "[INFO] Finished writing histograms.\n";
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
    Hist->Draw("HIST E1");
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
  //covarianceOsc *osc = nullptr;
  auto OscCovFile = GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
  auto OscCovName = GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");

  covarianceOsc* osc = new covarianceOsc(OscCovFile, OscCovName);

  auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});
  osc->setParameters(OscPars);
  
 // std::string hists_file (OutFileName + "event_histograms").c_str();
  // ####################################################################################
  // Create samplePDFFD objects

  std::vector<samplePDFFDBase *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

      // Create canvas
  auto gc1 = std::unique_ptr<TCanvas>(new TCanvas("gc1", "gc1", 800, 600));
  gStyle->SetOptStat(false);

  // Use only one multipage PDF for now
  std::string PdfFileName = "GenericBinTest.pdf";
  gc1->Print((PdfFileName + "[").c_str());  // open multipage PDF

  std::vector<TH1D*> DUNEHists;       // for PDFs
std::vector<TH1*> DUNEHistsall;     // for ROOT file

for (auto Sample : DUNEPdfs) {
    osc->setParameters(OscPars);

    // --- 1D histogram ---
    TH1D *nominalHist = (TH1D*)Sample->get1DHist()->Clone((Sample->GetTitle() + "_nominal").c_str());
    nominalHist->Sumw2();
    nominalHist->SetEntries(nominalHist->Integral());
    Sample->addData(nominalHist);
    Sample->reweight();

    if (xsec) {
        double nominal = xsec->getNominal(0);
        xsec->setParCurrProp(0, nominal);
        Sample->reweight();
    }

    TH1D *plotHist = (TH1D*)Sample->get1DHist()->Clone((Sample->GetTitle() + "_draw").c_str());
    plotHist->Sumw2();
    for (int i = 1; i <= plotHist->GetNbinsX(); i++)
        plotHist->SetBinError(i, std::sqrt(plotHist->GetBinContent(i)));
    plotHist->Scale(1, "WIDTH");

    DUNEHists.push_back(plotHist);     // for PDF
    DUNEHistsall.push_back(plotHist);  // for ROOT

    // Draw 1D histograms
    if (Sample->GetNDim() == 1) {
        plotHist->Draw("HIST E1");
        gc1->Print(PdfFileName.c_str());
    }

    // --- Generic binning ---
    if (Sample->generic_binning.GetNDimensions() > 0) {
        auto myhist = GetGenericBinningTH1(*Sample, Sample->GetTitle() + "_generic1D");
        myhist->Scale(1, "WIDTH");
        myhist->Draw();
        gc1->Print(PdfFileName.c_str());
        DUNEHistsall.push_back(myhist.release());  // transfer ownership

        if (Sample->generic_binning.GetNDimensions() == 2) {
            auto myhist2 = GetGenericBinningTH2(*Sample, Sample->GetTitle() + "_generic2D");
            myhist2->Draw("COLZ");
            gc1->Print(PdfFileName.c_str());
            DUNEHistsall.push_back(myhist2.release());

            auto slices1D = GetGenericBinningTH1Slices(*Sample, 0, Sample->GetTitle() + "_slice1D");
            int sliceIdx = 0;
            for (auto &slice : slices1D) {
                slice->SetName((Sample->GetTitle() + "_slice1D_" + std::to_string(sliceIdx++)).c_str());
                slice->Draw("COLZ");
                gc1->Print(PdfFileName.c_str());
                DUNEHistsall.push_back(slice.release());
            }
        }

        if (Sample->generic_binning.GetNDimensions() == 3) {
            auto slices2D = GetGenericBinningTH2Slices(*Sample, {0,1}, Sample->GetTitle() + "_slice2D");
            int sliceIdx = 0;
            for (auto &slice : slices2D) {
                slice->SetName((Sample->GetTitle() + "_slice2D_" + std::to_string(sliceIdx++)).c_str());
                slice->Draw("COLZ");
                gc1->Print(PdfFileName.c_str());
                DUNEHistsall.push_back(slice.release());
            }
        }
    }

    MACH3LOG_INFO("Event rate for {} : {:.2f}", Sample->GetTitle(), Sample->get1DHist()->Integral());
    MACH3LOG_INFO("LLH for {} : {:.10f}", Sample->GetTitle(), Sample->GetLikelihood());
}


  // Close the PDF properly
  gc1->Print((PdfFileName + "]").c_str());

  // Now write histograms
  WriteHistogramsToFile(OutFileName, DUNEHistsall);
  Write1DHistogramsToPdf(OutFileName, DUNEHists);

}

// #include <chrono>
// #include <iomanip>
// #include <iostream>
// #include <vector>
// #include <memory>

// #include <TCanvas.h>
// #include <TColor.h>
// #include <TH1D.h>
// #include <THStack.h>
// #include <TLegend.h>
// #include <TMath.h>
// #include <TRint.h>
// #include <TStyle.h>
// #include <TFile.h>
// #include <TROOT.h>
// #include <TSystem.h>

// #include "samplePDF/GenericBinningTools.h"
// #include "samplePDFDUNE/MaCh3DUNEFactory.h"
// #include "mcmc/mcmc.h"

// // ----------------------------------------------------------------------
// // Utility: Write TH1Ds to a ROOT file safely
// // ----------------------------------------------------------------------
// void Write1DHistogramsToFile(const std::string &OutFileName,
//                              const std::vector<TH1D *> &Histograms) {
//   std::cout << "\n[INFO] Writing histograms to ROOT file: " << OutFileName << std::endl;

//   TFile OutputFile(OutFileName.c_str(), "RECREATE");
//   if (!OutputFile.IsOpen()) {
//     std::cerr << "[ERROR] Could not open " << OutFileName << " for writing.\n";
//     return;
//   }

//   for (auto *Hist : Histograms) {
//     if (Hist) Hist->Write();
//   }

//   OutputFile.Write();   // flush
//   OutputFile.Close();   // close
//   std::cout << "[INFO] ROOT file successfully written.\n";
// }

// // ----------------------------------------------------------------------
// // Utility: Write TH1Ds to a multi-page PDF file
// // ----------------------------------------------------------------------
// void Write1DHistogramsToPdf(const std::string &BaseFileName,
//                             const std::vector<TH1D *> &Histograms) {
//   std::string PdfFileName = BaseFileName;
//   auto pos = PdfFileName.find(".root");
//   if (pos != std::string::npos) PdfFileName.replace(pos, 5, ".pdf");
//   else PdfFileName += ".pdf";

//   std::cout << "\n[INFO] Writing histograms to PDF: " << PdfFileName << std::endl;

//   auto c1 = std::make_unique<TCanvas>("c1", "c1", 800, 600);
//   c1->Print((PdfFileName + "[").c_str());  // open PDF
//   for (auto *Hist : Histograms) {
//     if (!Hist) continue;
//     c1->cd();
//     Hist->Draw("HIST E");
//     c1->Update(); // ensure rendering
//     c1->Print(PdfFileName.c_str());
//   }
//   c1->Print((PdfFileName + "]").c_str());  // close PDF

//   std::cout << "[INFO] PDF successfully written.\n";
// }

// // ----------------------------------------------------------------------
// // Main
// // ----------------------------------------------------------------------
// int main(int argc, char *argv[]) {
//   if (argc == 1) {
//     std::cout << "Usage: bin/EventRatesDUNEBeam config.cfg" << std::endl;
//     return 1;
//   }

//   // ------------------------------------------------------------------
//   // Setup and configuration
//   // ------------------------------------------------------------------
//   auto fitMan = std::make_unique<manager>(argv[1]);

//   std::string OutFileName = GetFromManager<std::string>(
//       fitMan->raw()["General"]["OutputFile"], "EventRates.root");

//   std::cout << "[INFO] Output ROOT file will be: " << OutFileName << std::endl;

//   // Derive corresponding PDF file name
//   std::string PrismFileName = OutFileName;
//   size_t pos = PrismFileName.find(".root");
//   if (pos != std::string::npos)
//       PrismFileName.replace(pos, 5, "_prism.pdf");
//   else
//       PrismFileName += "_prism.pdf";

//   std::cout << "[INFO] Prism PDF file will be: " << PrismFileName << std::endl;

//   // ------------------------------------------------------------------
//   // Initialize covariance objects
//   // ------------------------------------------------------------------
//   //covarianceXsec *xsec = new covarianceXsec();  // ✅ ensure initialized!
//   covarianceXsec *xsec = nullptr;
//   auto OscCovFile = GetFromManager<std::vector<std::string>>(
//       fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
//   auto OscCovName = GetFromManager<std::string>(
//       fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
//   auto *osc = new covarianceOsc(OscCovFile, OscCovName);

//   auto OscPars = GetFromManager<std::vector<double>>(
//       fitMan->raw()["General"]["OscillationParameters"], {});
//   osc->setParameters(OscPars);

//   // ------------------------------------------------------------------
//   // Create DUNE PDFs
//   // ------------------------------------------------------------------
//   std::vector<samplePDFFDBase *> DUNEPdfs;
//   MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

//   // Canvas for plots
//   auto gc1 = std::make_unique<TCanvas>("gc1", "gc1", 800, 600);
//   gStyle->SetOptStat(false);

//   gc1->Print((PrismFileName + "[").c_str());  // open multi-page PDF

//   std::vector<TH1D *> DUNEHists;

//   // ------------------------------------------------------------------
//   // Loop over samples
//   // ------------------------------------------------------------------
//   for (auto *Sample : DUNEPdfs) {
//     osc->setParameters(OscPars);
//     xsec->setParameters();
//     Sample->reweight();

//     // Make Asimov histogram
//     TH1D *Asimov_1D = (TH1D *)Sample->get1DHist()->Clone(
//         (Sample->GetTitle() + "_asimov").c_str());
//     Asimov_1D->SetEntries(Asimov_1D->Integral());
//     Sample->addData(Asimov_1D);

//     std::cout << "\n[INFO] Sample: " << Sample->GetTitle()
//               << " | Integral: " << Asimov_1D->Integral()
//               << " | Entries: " << Asimov_1D->GetEntries() << std::endl;

//     // Reweight & test xsec parameters
//     xsec->setParameters();
//     double nominal = xsec->getNominal(0);
//     double error = xsec->getDiagonalError(0);
//     xsec->setParCurrProp(0, nominal);

//     std::cout << "  Nominal parameter = " << nominal
//               << " ± " << error
//               << " (current: " << xsec->getParProp(0) << ")\n";

//     Sample->reweight();

//     // Draw and save
//     if (Sample->GetNDim() == 1) {
//       auto *h1 = (TH1D *)Sample->get1DHist()->Clone(
//           (Sample->GetTitle() + "_draw").c_str());
//       h1->Scale(1, "WIDTH");
//       h1->Sumw2();

//       h1->Draw("HIST E1");
//       gc1->Update();
//       gc1->Print(PrismFileName.c_str());
//       DUNEHists.push_back(h1);
//     }

//     DUNEHists.push_back(Sample->get1DHist());
//   }

//   // ------------------------------------------------------------------
//   // Close the multi-page PDF and write ROOT + summary PDF
//   // ------------------------------------------------------------------
//   gc1->Print((PrismFileName + "]").c_str());  // close multi-page PDF
//   gc1->Close();

//   Write1DHistogramsToFile(OutFileName, DUNEHists);
//   Write1DHistogramsToPdf(OutFileName, DUNEHists);

//   // ------------------------------------------------------------------
//   // Final ROOT cleanup
//   // ------------------------------------------------------------------
//   gROOT->GetListOfCanvases()->Delete();
//   gSystem->ProcessEvents();

//   std::cout << "\n✅ All files written successfully:\n"
//             << "   - " << OutFileName << "\n"
//             << "   - " << PrismFileName << "\n"
//             << "   - " << OutFileName.substr(0, OutFileName.find(".root")) + ".pdf" << "\n";

//   return 0;
// }
