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


// Returns all bin contents unfiltered (used to build mask and thrown predictions)
std::vector<double> GetObservablebinningUnfiltered(TH1* h, int ndim) {
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

  // -------------------------------------------------------
  // Fix all parameters
  // -------------------------------------------------------
  const int nPars = xsec->GetNParameters();
  for (int i = 0; i < nPars; ++i) {
    if (!xsec->IsParameterFixed(i)) {
      xsec->ToggleFixParameter(i);
    }
  }

  // -------------------------------------------------------
  // Unfix only DetSys parameters
  // -------------------------------------------------------
  auto xsecCovFiles = FitManager->raw()["General"]["Systematics"]["XsecCovFile"].as<std::vector<std::string>>();
  int nDetSys = 0;

  for (auto& covFile : xsecCovFiles) {
    YAML::Node covConfig = YAML::LoadFile(covFile);
    if (!covConfig["Systematics"]) continue;

    for (const auto& syst : covConfig["Systematics"]) {
      const auto& s = syst["Systematic"];
      if (!s["ParameterGroup"]) continue;
      if (s["ParameterGroup"].as<std::string>() != "DetSys") continue;

      std::string parName = s["Names"]["ParameterName"].as<std::string>();
      int idx = xsec->GetParIndex(parName);

      if (idx < 0) {
        MACH3LOG_WARN("DetSys parameter {} not found in handler, skipping", parName);
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
  if (nDetSys == 0) {
    MACH3LOG_ERROR("No DetSys parameters found! Check your XsecCovFile yamls and ParameterGroup labels.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Set all parameters to prior values before getting nominal
  for (int i = 0; i < nPars; ++i) {
    xsec->SetPar(i, xsec->GetParInit(i));
  }

  // -------------------------------------------------------
  // Get nominal prediction (unfiltered) and build mask
  // -------------------------------------------------------
  // The mask (validBins) is a list of flat indices into the unfiltered
  // concatenated bin vector. It is built once from the nominal and applied
  // identically to every throw, guaranteeing index alignment.

  std::vector<double> nominalFull;   // all bins, unfiltered
  std::vector<int>    sample_nbins_full;  // unfiltered bin counts per sample
  std::vector<int>    sample_nbins;       // masked bin counts per sample

  for (auto& pdf : DUNEPdfs) {
    pdf->Reweight();
    TH1* h = static_cast<TH1*>(pdf->GetMCHist(pdf->GetNDim())->Clone());
    std::vector<double> bins = GetObservablebinningUnfiltered(h, pdf->GetNDim());
    sample_nbins_full.push_back(bins.size());
    nominalFull.insert(nominalFull.end(), bins.begin(), bins.end());
    delete h;
  }

  // Build mask: keep only bins with non-zero nominal content
  std::vector<int> validBins;
  for (int i = 0; i < (int)nominalFull.size(); ++i) {
    if (nominalFull[i] > 0) validBins.push_back(i);
  }

  MACH3LOG_INFO("Total bins (unfiltered): {}", (int)nominalFull.size());
  MACH3LOG_INFO("Total bins (after masking empty bins): {}", (int)validBins.size());
  MACH3LOG_INFO("Bins removed (zero nominal content): {}", (int)nominalFull.size() - (int)validBins.size());

  // Lambda to apply mask to any full bin vector
  auto ApplyMask = [&](const std::vector<double>& full) {
    std::vector<double> masked;
    masked.reserve(validBins.size());
    for (int idx : validBins) masked.push_back(full[idx]);
    return masked;
  };

  // Masked nominal
  std::vector<double> nominal = ApplyMask(nominalFull);
  int totalBins = nominal.size();

  // Compute masked bin counts per sample (for boundary drawing later)
  {
    int flatIdx = 0;
    for (int s = 0; s < (int)DUNEPdfs.size(); ++s) {
      int count = 0;
      int start = flatIdx;
      int end   = flatIdx + sample_nbins_full[s];
      for (int idx : validBins)
        if (idx >= start && idx < end) ++count;
      sample_nbins.push_back(count);
      flatIdx += sample_nbins_full[s];
    }
  }

  // -------------------------------------------------------
  // Initialise covariance accumulator
  // -------------------------------------------------------
  TMatrixDSym CovMatrix(totalBins);
  CovMatrix.Zero();

  // -------------------------------------------------------
  // Throw loop
  // -------------------------------------------------------
  std::vector<TH1D*> histogram_per_bin;
  for (int i = 0; i < totalBins; ++i) {
    std::string hname = "bin_" + std::to_string(i);
    histogram_per_bin.push_back(new TH1D(hname.c_str(), hname.c_str(), 100, +1, -1));
  }

  for (int iThrow = 0; iThrow < nThrows; ++iThrow) {

    xsec->ThrowParameters();

    // Get full unfiltered prediction for this throw
    std::vector<double> thrownFull;
    thrownFull.reserve(nominalFull.size());

    for (auto& pdf : DUNEPdfs) {
      pdf->Reweight();
      TH1* h = static_cast<TH1*>(pdf->GetMCHist(pdf->GetNDim())->Clone());
      std::vector<double> bins = GetObservablebinningUnfiltered(h, pdf->GetNDim());
      thrownFull.insert(thrownFull.end(), bins.begin(), bins.end());
      delete h;
    }

    // Apply the same mask as the nominal — guarantees index alignment
    std::vector<double> thrown_pred = ApplyMask(thrownFull);

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

    if (iThrow % 100 == 0) {
      MACH3LOG_INFO("Throw {}/{}", iThrow, nThrows);
    }
  }

  // Normalise
  double norm = 1.0 / (double)nThrows;
  for (int i = 0; i < totalBins; ++i)
    for (int j = 0; j < totalBins; ++j)
      CovMatrix(i, j) *= norm;

  // -------------------------------------------------------
  // Fractional covariance
  // -------------------------------------------------------
  TMatrixDSym FracCovMatrix(totalBins);
  for (int i = 0; i < totalBins; ++i) {
    for (int j = 0; j < totalBins; ++j) {
      if (nominal[i] > 0 && nominal[j] > 0)
        FracCovMatrix(i, j) = CovMatrix(i, j) / (nominal[i] * nominal[j]);
      else
        FracCovMatrix(i, j) = 0.0;
    }
  }

  // -------------------------------------------------------
  // Correlation matrix
  // -------------------------------------------------------
  TMatrixDSym CorrMatrix(totalBins);
  for (int i = 0; i < totalBins; ++i) {
    for (int j = 0; j < totalBins; ++j) {
      double denom = std::sqrt(CovMatrix(i,i) * CovMatrix(j,j));
      CorrMatrix(i, j) = denom > 0 ? CovMatrix(i, j) / denom : 0.0;
    }
  }

  // Sanity check: all diagonal entries should be 1.0
  MACH3LOG_INFO("Checking correlation matrix diagonal (all should be 1.0):");
  int nBadDiag = 0;
  for (int i = 0; i < totalBins; ++i) {
    if (std::abs(CorrMatrix(i,i) - 1.0) > 1e-6) {
      MACH3LOG_WARN("  bin {} diagonal = {:.6f} (expected 1.0)", i, CorrMatrix(i,i));
      ++nBadDiag;
    }
  }
  if (nBadDiag == 0)
    MACH3LOG_INFO("  All diagonal entries are 1.0.");

  // Print fractional uncertainty per bin
  MACH3LOG_INFO("Fractional uncertainty (sqrt of frac cov diagonal) per bin:");
  for (int i = 0; i < totalBins; ++i) {
    if (FracCovMatrix(i,i) > 0)
      MACH3LOG_INFO("  bin {}: {:.2f}%", i, std::sqrt(FracCovMatrix(i,i)) * 100.0);
  }

  // -------------------------------------------------------
  // Save output
  // -------------------------------------------------------
  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  CovMatrix.Write("NDCovMatrix");
  FracCovMatrix.Write("NDFracCovMatrix");
  CorrMatrix.Write("NDCorrMatrix");

  TH1D hNominal("Nominal", "Nominal", totalBins, 0, totalBins);
  for (int i = 0; i < totalBins; ++i)
    hNominal.SetBinContent(i+1, nominal[i]);
  hNominal.Write();

  // Save the valid bin index map so consumers know which original bins are included
  TH1I hValidBinMap("ValidBinMap", "Flat index of each masked bin in the unfiltered vector",
                    totalBins, 0, totalBins);
  for (int i = 0; i < totalBins; ++i)
    hValidBinMap.SetBinContent(i+1, validBins[i]);
  hValidBinMap.Write();

  TH1D hBinMap("BinMap", "Sample bin ranges", DUNEPdfs.size(), 0, DUNEPdfs.size());
  int binOffset = 0;
  for (unsigned s = 0; s < DUNEPdfs.size(); ++s) {
    MACH3LOG_INFO("Sample {}: {} masked bins (bins {} to {})",
                  DUNEPdfs[s]->GetTitle(), sample_nbins[s],
                  binOffset, binOffset + sample_nbins[s] - 1);
    hBinMap.GetXaxis()->SetBinLabel(s+1, DUNEPdfs[s]->GetTitle().c_str());
    hBinMap.SetBinContent(s+1, sample_nbins[s]);
    binOffset += sample_nbins[s];
  }
  hBinMap.Write();

  // -------------------------------------------------------
  // Plot matrices to PDF
  // -------------------------------------------------------
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
  hCorr->SetMaximum( 1.0);

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
      line.DrawLine(offset, 0,         offset,    totalBins);
      line.DrawLine(0,      offset,    totalBins, offset);
    }
  };

  c->Print((pdfName + "[").c_str());

  hCov->SetMinimum(-10000.0);
  hCov->SetMaximum( 10000.0);
  hCov->Draw("COLZ");
  DrawLines();
  c->Print(pdfName.c_str());

  hFrac->SetMinimum(0.0);
  hFrac->SetMaximum(1e-3);
  hFrac->Draw("COLZ");
  DrawLines();
  c->Print(pdfName.c_str());

  hCorr->Draw("COLZ");
  DrawLines();
  c->Print(pdfName.c_str());

  c->Print((pdfName + "]").c_str());
  MACH3LOG_INFO("Matrix plots written to {}", pdfName);

  // Save per-bin throw distributions
  TDirectory* throwDir = OutputFile->mkdir("ThrowDistributions");
  throwDir->cd();
  for (int i = 0; i < totalBins; ++i)
    histogram_per_bin[i]->Write();
  OutputFile->cd();

  // Cleanup
  for (auto h : histogram_per_bin) delete h;
  delete hCov; delete hFrac; delete hCorr; delete c;
  OutputFile->Close();

  MACH3LOG_INFO("Done. Covariance matrix written to {}", OutputFileName);

  return 0;
}