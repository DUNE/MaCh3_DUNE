#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"

#include "TFile.h"

#include "yaml-cpp/yaml.h"

#include "eigen/Dense"

#include <cmath>
#include <iostream>
#include <vector>

// #define DEBUG_MAKE_COV_MX

int main(int argc, char *argv[]) {

  if (argc <= 5) {
    std::cout << "[ERROR]: Runlike " << argv[0]
              << " <Mach3_Config.yml> <Parameter.yml> <NThrows> <output.root>"
              << std::endl;
  }

  auto conf_yml = M3OpenConfig(argv[1]);
  conf_yml["General"]["Systematics"] =
      YAML::Load(std::string("XsecCovFile:\n\t[") + argv[2] +
                 "]\nXsecCovName: xsec_cov\nXsecStepScale: 1.0");
  auto fitMan = std::make_unique<Manager>(conf_yml);

  auto nThrows = std::stol(argv[3]);
  auto OutputFileName = std::string(argv[4]);

  ParameterHandlerGeneric *xsec = nullptr;
  std::vector<SampleHandlerFD *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  if (DUNEPdfs.size() != 1) {
    throw std::runtime_error("only want a single samplehandler");
  }

  auto &pdf = *DUNEPdfs.front();

  auto get_pred = [&]() -> Eigen::VectorXd {
    pdf.Reweight();
    auto mcvals = pdf.GetMCArray();

#ifdef DEBUG_MAKE_COV_MX
    std::cout << "[ ";
    for (auto &v : mcvals) {
      std::cout << v << ", ";
    }
    std::cout << "]" << std::endl;
#endif

    return Eigen::Map<Eigen::VectorXd>(mcvals.data(), mcvals.size());
  };

  // Set all parameters to prior values before getting nominal
  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    xsec->SetPar(i, xsec->GetParInit(i));
  }
  Eigen::VectorXd cvpred = get_pred();
  Eigen::VectorXd mean = Eigen::VectorXd::Zero(cvpred.size());

  Eigen::MatrixXd CovMatrix =
      Eigen::MatrixXd::Zero(cvpred.size(), cvpred.size());

  for (int iThrow = 0; iThrow < nThrows; ++iThrow) {

    xsec->ThrowParameters();

    Eigen::VectorXd throw_pred = get_pred();

    mean += throw_pred;

    Eigen::VectorXd diff = (throw_pred - cvpred);

#ifdef DEBUG_MAKE_COV_MX
    std::cout << "diff: [ ";
    for (int i = 0; i < diff.rows(); ++i) {
      std::cout << diff(i, 0) << ", ";
    }
    std::cout << "]" << std::endl;
#endif

    CovMatrix += (diff * diff.transpose());

    if (iThrow % 100 == 0) {
      MACH3LOG_INFO("Throw {}/{}", iThrow, nThrows);
    }
  }

  mean *= 1.0 / double(nThrows);
  CovMatrix *= 1.0 / double(nThrows - 1);

  auto OutputFile =
      std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  TMatrixD cvpred_root(cvpred.size(), 1, cvpred.data());
  TMatrixD mean_root(mean.size(), 1, mean.data());
  TMatrixD CovMatrix_root(CovMatrix.rows(), CovMatrix.cols(), CovMatrix.data());

  cvpred_root.Write("NDCV");
  mean_root.Write("NDThrowsMean");
  CovMatrix_root.Write("NDCovMatrix");

  // write PDF of covmatrix (correlation matrix)

  OutputFile->Write();
  OutputFile->Close();
}
