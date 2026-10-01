#pragma once

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
