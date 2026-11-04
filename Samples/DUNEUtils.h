#pragma once

#include "Samples/BinningHandler.h"

#include <iomanip>
#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <openssl/evp.h>

template<class Matrix>
inline void WriteEigenMatrixToFile(const std::string& filename, const Matrix& matrix){
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  if(out.is_open()) {
    typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
    out.write(reinterpret_cast<char*>(&rows), sizeof(typename Matrix::Index));
    out.write(reinterpret_cast<char*>(&cols), sizeof(typename Matrix::Index));
    out.write(reinterpret_cast<const char*>(matrix.data()), rows*cols*static_cast<typename Matrix::Index>(sizeof(typename Matrix::Scalar)) );
    out.close();
  }
  else {
    MACH3LOG_ERROR("Can not write to file:{}",filename);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

template<class Matrix>
inline void ReadEigenMatrixFromFile(const std::string& filename, Matrix& matrix){
  std::ifstream in(filename, std::ios::in | std::ios::binary);
  if (in.is_open()) {
    typename Matrix::Index rows=0, cols=0;
    in.read(reinterpret_cast<char*>(&rows),sizeof(typename Matrix::Index));
    in.read(reinterpret_cast<char*>(&cols),sizeof(typename Matrix::Index));
    matrix.resize(rows, cols);
    in.read(reinterpret_cast<char*>(matrix.data()), rows*cols*static_cast<typename Matrix::Index>(sizeof(typename Matrix::Scalar)) );
    in.close();
  }
  else {
    MACH3LOG_ERROR("Can not open binary matrix file:{}",filename);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// Compute MD5 digest of buffer 'data'. Returns true on success and fills 'digest' (16 bytes).
inline void ComputeMD5(const unsigned char* data, size_t len, unsigned char digest[EVP_MAX_MD_SIZE], unsigned int &digest_len) {
  EVP_MD_CTX *ctx = nullptr;
  
  ctx = EVP_MD_CTX_new();
  if (!ctx) {
    std::cerr << "EVP_MD_CTX_new failed\n";
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  // Use EVP_md5() — still available under OpenSSL 3.0's default provider.
  if (1 != EVP_DigestInit_ex(ctx, EVP_md5(), nullptr)) {
    std::cerr << "EVP_DigestInit_ex failed\n";
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  if (len > 0) {
    if (1 != EVP_DigestUpdate(ctx, data, len)) {
      std::cerr << "EVP_DigestUpdate failed\n";
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  
  if (1 != EVP_DigestFinal_ex(ctx, digest, &digest_len)) {
    std::cerr << "EVP_DigestFinal_ex failed\n";
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  EVP_MD_CTX_free(ctx);
}

// Read entire file into vector<char>, returns false on error
inline bool ReadFileIntoVector(const std::string &path, std::vector<unsigned char> &out) {
  std::ifstream ifs(path, std::ios::binary | std::ios::ate);
  if (!ifs) return false;
  std::streamsize size = ifs.tellg();
  ifs.seekg(0, std::ios::beg);
  out.resize(static_cast<size_t>(size));
  if (!ifs.read(reinterpret_cast<char*>(out.data()), size)) return false;
  return true;
}

// Convert binary digest to lowercase hex string
inline std::string ToHex(const unsigned char* data, size_t len) {
  std::ostringstream oss;
  oss << std::hex << std::setfill('0');
  for (size_t i = 0; i < len; ++i) {
    oss << std::setw(2) << static_cast<int>(data[i]);
  }
  return oss.str();
}

inline bool CompareMD5Sum(std::string KnownSum, std::string FilePathToCheck) {
  MACH3LOG_INFO("Determining MD5 Sum for File:{}",FilePathToCheck);
  
  std::vector<unsigned char> Data;
  if (!ReadFileIntoVector(FilePathToCheck,Data)) {
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  unsigned char Digest[EVP_MAX_MD_SIZE];
  unsigned int Digest_len = 0;
  ComputeMD5(Data.data(), Data.size(), Digest, Digest_len);
  std::string CalculatedMD5Sum = ToHex(Digest, Digest_len);

  MACH3LOG_INFO("Comparing to known MD5 Sum: {}",KnownSum);
  if (KnownSum != CalculatedMD5Sum) {
    MACH3LOG_ERROR("MD5 Sums are different! Calculated MD5: {}",CalculatedMD5Sum);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  MACH3LOG_INFO("MD5 Sums are identical!");
  return true;
}

inline std::vector<std::unique_ptr<TH2D>> Make2DSlices(const TH1* hist, const BinningHandler& binning, int iSample, std::vector<std::string> BinningVariables, int sliceAxis = 0) {
  if (!hist)
    throw std::invalid_argument("Null TH1 pointer");
  
  if (binning.GetNDim(iSample) != 3)
    throw std::runtime_error("Expected 3D binning");
  
  if (!binning.IsUniform(iSample))
    throw std::runtime_error("Requires Cartesian (Uniform) binning");
  
  if (sliceAxis < 0 || sliceAxis >= 3)
    throw std::invalid_argument("Invalid slice axis");
  
  if (hist->GetNbinsX() != binning.GetNBins(iSample))
    throw std::runtime_error("TH1 bin count mismatch");
  
  // Determine which axes to display.
  std::array<int, 2> axes;
  int a = 0;
  
  for (int d = 0; d < 3; ++d) {
    if (d != sliceAxis)
      axes[a++] = d;
  }
  
  const int xAxis = axes[0];
  const int yAxis = axes[1];
  
  const auto xEdges = binning.GetBinEdges(iSample, xAxis);
  const auto yEdges = binning.GetBinEdges(iSample, yAxis);
  const auto sEdges = binning.GetBinEdges(iSample, sliceAxis);
  
  const int nx = binning.GetNAxisBins(iSample, xAxis);
  const int ny = binning.GetNAxisBins(iSample, yAxis);
  const int ns = binning.GetNAxisBins(iSample, sliceAxis);
  
  std::vector<std::unique_ptr<TH2D>> slices;
  slices.reserve(ns);
  
  for (int s = 0; s < ns; ++s) {
    
    const TString name = Form("%s_slice%d_bin%d", hist->GetName(), sliceAxis, s);
    const TString title = Form("%s: %.3g - %.3g;%s;%s", BinningVariables[sliceAxis].c_str(), sEdges[s], sEdges[s+1], BinningVariables[xAxis].c_str(), BinningVariables[yAxis].c_str());

    auto h2 = std::make_unique<TH2D>(name, title, nx, xEdges.data(), ny, yEdges.data());

    h2->SetDirectory(nullptr);
    
    for (int ix = 0; ix < nx; ++ix) {
      for (int iy = 0; iy < ny; ++iy) {
	
	std::vector<int> indices(3);
	
	indices[sliceAxis] = s;
	indices[xAxis] = ix;
	indices[yAxis] = iy;
	
	// MaCh3 returns zero-based sample-local bin.
	const int flatBin = binning.GetBinSafe(iSample, indices);
	
	// ROOT TH1 bins start at 1.
	const int rootBin = flatBin + 1;
	
	h2->SetBinContent(ix + 1, iy + 1, hist->GetBinContent(rootBin));
      }
    }
    
    slices.push_back(std::move(h2));
  }
  
  return slices;
}

std::unique_ptr<TCanvas> Plot2DSlices(const std::vector<std::unique_ptr<TH2D>>& slices, const std::string& canvasName = "Slices", bool LogX = false, bool LogY = false, int nCols = 5) {
  const int nSlices = static_cast<int>(slices.size());
  
  if (nSlices == 0)
    return nullptr;
  
  nCols = std::max(1, std::min(nCols, nSlices));
  const int nRows = (nSlices + nCols - 1) / nCols;
  
  auto canvas = std::make_unique<TCanvas>(canvasName.c_str(),canvasName.c_str(),350 * nCols,320 * nRows);
  canvas->Divide(nCols, nRows, 0.005, 0.005);
  
  // Find global Z range across all slices.
  double zMin = slices[0]->GetMinimum();
  double zMax = slices[0]->GetMaximum();
  
  for (const auto& hist : slices) {
    zMin = std::min(zMin, hist->GetMinimum());
    zMax = std::max(zMax, hist->GetMaximum());
  }
  
  gStyle->SetOptStat(0);
  
  for (int i = 0; i < nSlices; ++i) {
    canvas->cd(i + 1);

    gPad->SetLogx(LogX);
    gPad->SetLogy(LogY);
    
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.12);
    
    auto* hist = slices[i].get();
    
    hist->SetMinimum(zMin);
    hist->SetMaximum(zMax);
    
    hist->GetXaxis()->SetTitleSize(0.045);
    hist->GetYaxis()->SetTitleSize(0.045);
    hist->GetZaxis()->SetTitleSize(0.045);
    
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.035);
    hist->GetZaxis()->SetLabelSize(0.035);
    
    hist->Draw("COLZ");
  }
  
  canvas->cd();
  canvas->Update();
  
  return canvas;
}
