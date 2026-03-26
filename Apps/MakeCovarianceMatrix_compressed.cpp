#include <iostream>
#include <vector>
#include <cmath>
#include <TH1D.h>
#include <TFile.h>
#include <TMatrixDSym.h>
#include "yaml-cpp/yaml.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLine.h>
#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"


std::vector<double> GetObservablebinning(TH1* h, int ndim) {
  std::vector<double> bins;
  if (ndim == 2) {
    TH2* h2 = static_cast<TH2*>(h);
    for (int ix = 1; ix <= h2->GetNbinsX(); ++ix)
      for (int iy = 1; iy <= h2->GetNbinsY(); ++iy)
        bins.push_back(h2->GetBinContent(ix, iy));
  } else {
    for (int i = 1; i <= h->GetNbinsX(); ++i)
      bins.push_back(h->GetBinContent(i));
  }
  return bins;
}



int main(int argc, char* argv[]) {

  const int NRGBs = 3;
  const int NCont = 255;
  double stops[NRGBs] = { 0.0, 0.5, 1.0 };
  double red[NRGBs]   = { 0.0, 1.0, 1.0 };
  double green[NRGBs] = { 0.0, 1.0, 0.0 };
  double blue[NRGBs]  = { 1.0, 1.0, 0.0 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);


  auto FitManager = MaCh3ManagerFactory(argc, argv);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();
  int nThrows = GetFromManager<int>(FitManager->raw()["General"]["CovarianceMatrixParams"]["NThrows"], 1000);

  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);


  // Fix all parameters
  const int nPars = xsec->GetNParameters();
  for (int i = 0; i < nPars; ++i) {
    if (!xsec->IsParameterFixed(i)) {
      xsec->ToggleFixParameter(i);
    }
  }

  // Free only the detector systematic parameters
  auto xsecCovFiles = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();
  int nDetSys = 0;

  for (auto& covFile : xsecCovFiles) {
    YAML::Node covConfig = YAML::LoadFile(covFile);
    if (!covConfig["Systematics"]) continue;

    for (const auto& syst : covConfig["Systematics"]) {
      const auto& s = syst["Systematic"];
      if (!s["ParameterGroup"]) continue;
      if (s["ParameterGroup"].as<std::string>() != "Xsec") continue;

      std::string parName = s["Names"]["ParameterName"].as<std::string>();
      int idx = xsec->GetParIndex(parName);

      if (idx < 0) {
        MACH3LOG_WARN("Xsec parameter {} not found in handler, skipping", parName);
        continue;
      }

      if (xsec->IsParameterFixed(idx)) {
        xsec->ToggleFixParameter(idx);
        MACH3LOG_INFO("Throwing detector systematic: {} (idx {})", parName, idx);
        nDetSys++;
      } else {
        MACH3LOG_WARN("DetSys parameter {} was already unfixed - not toggling", parName);
      }
    }
  }

  MACH3LOG_INFO("Total DetSys parameters being thrown: {}", nDetSys);

  // Set all parameters to prior values before getting nominal
  for (int i = 0; i < nPars; ++i) {
    xsec->SetPar(i, xsec->GetParInit(i));
  }


  // Get nominal prediction — only keep bins with events
  std::vector<double> nominal;
  std::vector<int> sample_nbins;
  // Maps each  bin index (bins with events in  to the sample index, local bin index within that sample)
  std::vector<std::pair<int,int>> activeBinOrigin;

  for (unsigned s = 0; s < DUNEPdfs.size(); ++s) {
    auto& pdf = DUNEPdfs[s];
    pdf->Reweight();
    TH1* h = static_cast<TH1*>(pdf->GetMCHist(pdf->GetNDim())->Clone());
    std::vector<double> bins = GetObservablebinning(h, pdf->GetNDim());

    int bin_haseventsin = 0;
    for (int ib = 0; ib < (int)bins.size(); ++ib) {
      if (bins[ib] > 0.0) {
        nominal.push_back(bins[ib]);
        activeBinOrigin.push_back({(int)s, ib});
        ++bin_haseventsin;
      }
    }
    sample_nbins.push_back(bin_haseventsin);
    MACH3LOG_INFO("Sample {}: {}/{} bins has non-zero numer of events in",
                  pdf->GetTitle(), bin_haseventsin, (int)bins.size());
    delete h;
  }

  int totalBins = nominal.size();
  MACH3LOG_INFO("Total bins containing events across all samples: {}", totalBins);


  // Initialise covariance accumulator
  TMatrixDSym CovMatrix(totalBins);
  CovMatrix.Zero();


  // Initialise per-bin throw histograms
  std::vector<TH1D*> histogram_per_bin;
  for (int i = 0; i < totalBins; ++i) {
    std::string hname = "bin_" + std::to_string(i);
    double nom = nominal[i];
    double range = std::max(nom * 0.5, 1.0);
    histogram_per_bin.push_back(
      new TH1D(hname.c_str(), hname.c_str(), 100, nom - range, nom + range));
  }


  // Throw loop
  for (int iThrow = 0; iThrow < nThrows; ++iThrow) {

    xsec->ThrowParameters();
    xsec->PrintFunctionalParams();

    // Collect all bins from all samples for this throw, indexed by sample
    std::vector<std::vector<double>> allSampleBins(DUNEPdfs.size());
    for (unsigned s = 0; s < DUNEPdfs.size(); ++s) {
      auto& pdf = DUNEPdfs[s];
      pdf->Reweight();
      TH1* h = static_cast<TH1*>(pdf->GetMCHist(pdf->GetNDim())->Clone());
      allSampleBins[s] = GetObservablebinning(h, pdf->GetNDim());
      delete h;
    }

    // Build thrown prediction using the same active bins as the nominal
    std::vector<double> thrown_pred;
    thrown_pred.reserve(totalBins);
    for (auto& [si, ib] : activeBinOrigin) {
      thrown_pred.push_back(allSampleBins[si][ib]);
    }

    // Fill per-bin histograms
    for (int i = 0; i < totalBins; ++i) {
      histogram_per_bin[i]->Fill(thrown_pred[i]);
    }

    // Accumulate covariance
    for (int i = 0; i < totalBins; ++i) {
      double delta_i = thrown_pred[i] - nominal[i];
      for (int j = i; j < totalBins; ++j) {
        double delta_j = thrown_pred[j] - nominal[j];
        CovMatrix(i, j) += delta_i * delta_j;
        if (i != j) CovMatrix(j, i) += delta_i * delta_j;
      }
    }

   
    const int printpredictions = 20; 
    std::cout << "Bin | Nominal | Thrown | Delta\n";
    std::cout << "-------------------------------------------------\n";
    for (int i = 0; i < std::min(printpredictions, totalBins); ++i) {
        double delta = thrown_pred[i] - nominal[i];
        std::cout << std::setw(4)  << i             << " | "
                  << std::setw(10) << nominal[i]     << " | "
                  << std::setw(10) << thrown_pred[i] << " | "
                  << std::setw(10) << delta           << "\n";
    }

    if (iThrow % 100 == 0) {
      MACH3LOG_INFO("Throw {}/{}", iThrow, nThrows);
    }
  }

  


  // Normalise cov matrix
  double norm = 1.0 / (double)nThrows;
  for (int i = 0; i < totalBins; ++i)
    for (int j = 0; j < totalBins; ++j)
      CovMatrix(i, j) *= norm;
      


  // Fractional covariance
  TMatrixDSym FracCovMatrix(totalBins);
  for (int i = 0; i < totalBins; ++i) {
    for (int j = 0; j < totalBins; ++j) {
      if (nominal[i] > 0 && nominal[j] > 0)
        FracCovMatrix(i, j) = CovMatrix(i, j) / (nominal[i] * nominal[j]);
      else
        FracCovMatrix(i, j) = 0.0;
    }
  }

  // Correlation matrix — diagonal is explicitly set to 1.0
  TMatrixDSym CorrMatrix(totalBins);
  for (int i = 0; i < totalBins; ++i) {
    for (int j = 0; j < totalBins; ++j) {
      if (i == j) {
        CorrMatrix(i, j) = 1.0;
      } else {
        double denom = std::sqrt(CovMatrix(i,i) * CovMatrix(j,j));
        CorrMatrix(i, j) = denom > 0 ? CovMatrix(i, j) / denom : 0.0;
      }
    }
  }

  for (int i = 0; i < totalBins; ++i) {
    for (int j = 0; j < totalBins; ++j) {
      std::cout << "Cov matrix element (i,j) = " << CovMatrix(i,j) << std::endl;
    }
  }
  // Print fractional uncertainty per bin for sanity check
  MACH3LOG_INFO("Fractional uncertainty (sqrt of frac cov diagonal) per bin:");
  for (int i = 0; i < totalBins; ++i) {
    if (FracCovMatrix(i,i) > 0)
      MACH3LOG_INFO("  bin {}: {:.2f}%", i, std::sqrt(FracCovMatrix(i,i)) * 100.0);
  }


  // Write outputs
  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  CovMatrix.Write("NDCovMatrix");
  FracCovMatrix.Write("NDFracCovMatrix");
  CorrMatrix.Write("NDCorrMatrix");

  TH1D hNominal("Nominal", "Nominal", totalBins, 0, totalBins);
  for (int i = 0; i < totalBins; ++i)
    hNominal.SetBinContent(i+1, nominal[i]);
  hNominal.Write();

  TH1D hBinMap("BinMap", "Sample bin ranges", DUNEPdfs.size(), 0, DUNEPdfs.size());
  int binOffset = 0;
  for (unsigned s = 0; s < DUNEPdfs.size(); ++s) {
    MACH3LOG_INFO("Sample {}: bins {} to {}", DUNEPdfs[s]->GetTitle(), binOffset, binOffset + sample_nbins[s] - 1);
    hBinMap.GetXaxis()->SetBinLabel(s+1, DUNEPdfs[s]->GetTitle().c_str());
    hBinMap.SetBinContent(s+1, sample_nbins[s]);
    binOffset += sample_nbins[s];
  }
  hBinMap.Write();

  gStyle->SetOptStat(0);

  auto MatrixToHist = [&](const TMatrixDSym& mat, const char* name, const char* title) -> TH2D* {
    TH2D* h = new TH2D(name, title, totalBins, 0, totalBins, totalBins, 0, totalBins);
    for (int i = 0; i < totalBins; ++i)
      for (int j = 0; j < totalBins; ++j)
        h->SetBinContent(i+1, j+1, mat(i,j));
    return h;
  };

  auto DrawSampleBoundaries = [&](TH2D* h) {
    int offset = 0;
    for (unsigned s = 0; s < DUNEPdfs.size(); ++s) {
      int mid = offset + sample_nbins[s] / 2;
      h->GetXaxis()->SetBinLabel(mid + 1, DUNEPdfs[s]->GetTitle().c_str());
      h->GetYaxis()->SetBinLabel(mid + 1, DUNEPdfs[s]->GetTitle().c_str());
      offset += sample_nbins[s];
    }
  };

  TH2D* hCov  = MatrixToHist(CovMatrix,    "hCovMatrix",     "ND Covariance Matrix");
  TH2D* hFrac = MatrixToHist(FracCovMatrix, "hFracCovMatrix", "ND Fractional Covariance Matrix");
  TH2D* hCorr = MatrixToHist(CorrMatrix,    "hCorrMatrix",    "ND Correlation Matrix");

  DrawSampleBoundaries(hCov);
  DrawSampleBoundaries(hFrac);
  DrawSampleBoundaries(hCorr);

  hCorr->SetMinimum(-1.0);
  hCorr->SetMaximum(1.0);

  std::string pdfName = OutputFileName.substr(0, OutputFileName.find_last_of('.')) + "_matrices.pdf";
  TCanvas* c = new TCanvas("c", "c", 1200, 1000);
  c->SetRightMargin(0.15);

  auto DrawLines = [&]() {
    TLine line;
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.SetLineColor(kRed);
    int offset = 0;
    for (unsigned s = 0; s < DUNEPdfs.size() - 1; ++s) {
      offset += sample_nbins[s];
      line.DrawLine(offset, 0, offset, totalBins);
      line.DrawLine(0, offset, totalBins, offset);
    }
  };

  c->Print((pdfName + "[").c_str());

  
  double covMax = 0.0;
  double fracMax = 0.0;
  for (int i = 0; i < totalBins; ++i)
    for (int j = 0; j < totalBins; ++j)
      covMax = std::max(covMax, std::abs(CovMatrix(i, j)));
     

  hCov->SetMinimum(0);
  hCov->SetMaximum(200000);

  
  hFrac->SetMinimum(0.0);
  hFrac->SetMaximum(1e-2);

  hCov->Draw("COLZ");
  DrawLines();
  c->Print(pdfName.c_str());

  hFrac->Draw("COLZ");
  DrawLines();
  c->Print(pdfName.c_str());

  hCorr->Draw("COLZ");
  DrawLines();
  c->Print(pdfName.c_str());

  // Save per-bin throw distribution to see if throws look Gaussian
  TDirectory* throwDir = OutputFile->mkdir("ThrowDistributions");
  throwDir->cd();
  for (int i = 0; i < totalBins; ++i) {
    histogram_per_bin[i]->Write();
  }
  OutputFile->cd();

  
  for (auto h : histogram_per_bin) delete h;

  c->Print((pdfName + "]").c_str());


  delete hCov; delete hFrac; delete hCorr; delete c;
  OutputFile->Close();

  MACH3LOG_INFO("Covariance matrix written to {}", OutputFileName);

  return 0;
}